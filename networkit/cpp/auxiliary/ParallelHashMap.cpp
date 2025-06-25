#include <omp.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cuchar>
#include <iterator>
#include <thread>

#include <networkit/auxiliary/ParallelHashMap.hpp>

namespace Aux {

namespace {
// From 5.3.1 of "Concurrent Hash Tables: Fast and general(?)!" from Maier et
// al. (2016): h_c(x) := floor(h(x) * (c/U)) where c = capacity and U = key
// space = 64 (bit) in our case. Since capacity and U are both power of 2 we can
// change the formula to use bit-shifting.
uint64_t scaling_log(uint64_t const h, uint64_t const log_2_capacity) {
    return h >> (ParallelHashMap::ht_key_space - log_2_capacity);
}

bool ht_filled(size_t const occupancy, size_t const capacity) {
    uint64_t const alpha_elements_count = capacity >> 1;
    return occupancy > alpha_elements_count;
}

ParallelHashMap::HTAtomic128 *grow_hashtable(std::unique_ptr<ParallelHashMap::HTAtomic128> &source,
                                             std::unique_ptr<ParallelHashMap::HTAtomic128> &target,
                                             std::atomic_bool &request_growth, int const p_count,
                                             int const p_id) {
    request_growth.store(true);

#pragma omp single
    { target = std::make_unique<ParallelHashMap::HTAtomic128>(source->capacity() << 1); }
    // omp implicit barrier

    source->roam(*target, p_count, p_id);
#pragma omp barrier

#pragma omp single
    {
        request_growth.store(false);
        source = std::move(target);
    }
    // omp implicit barrier

    return source.get();
}

} // namespace

using HTAtomic128 = ParallelHashMap::HTAtomic128;
using HTSyncData = ParallelHashMap::HTSyncData;
using HTHandle = ParallelHashMap::HTHandle;

HTAtomic128::HTAtomic128() : m_cells(m_capacity, Cell{}) {}

HTAtomic128::HTAtomic128(size_t const capacity)
    : m_capacity(capacity), m_log_capacity(std::log2(capacity)), m_cells(capacity, Cell{}) {}

HTAtomic128::HTAtomic128(HTAtomic128 const &other)
    : m_capacity(other.m_capacity), m_log_capacity(other.m_log_capacity), m_cells(other.m_cells),
      m_global_occupancy(other.m_global_occupancy.load()) {}

HTAtomic128::HTAtomic128(HTAtomic128 &&other) noexcept : HTAtomic128() {
    swap(*this, other);
}

HTAtomic128 &HTAtomic128::operator=(HTAtomic128 other) {
    swap(*this, other);
    return *this;
}

bool HTAtomic128::insert(uint64_t const key, uint64_t const value) {
    Cell const desired = makeCell(key, value);
    Cell expected = invalidCell();
    uint64_t const hashed_key = hash64(key);
    size_t idx = scaling_log(hashed_key, m_log_capacity);

    while (true) {
        idx = fastMod(idx, m_capacity);
        assert(idx < m_capacity);

        Cell &cell = m_cells[idx];

#pragma omp atomic read
        expected.key = cell.key;

        if (desired.key == expected.key) {
            return false;
        }

#pragma omp atomic read
        expected.value = cell.value;

        // Data Race check on invalid and desired
        if (expected.key == ParallelHashMap::ht_invalid_key) {
            if (dwcas(cell, expected, desired)) {
                return true;
            }
            // DWCAS did not work.. thus reaching this line due to a race on
            // the cells index with another thread when trying to write. We
            // simply don't increment the idx variable which will try again
            // with the same idx and proceed with linear probing.
        } else {
            idx++;
        }
    }
}

uint64_t HTAtomic128::find(uint64_t const key) const {
    Cell actual = invalidCell();
    uint64_t const hashed_key = hash64(key);
    size_t idx = scaling_log(hashed_key, m_log_capacity);

    while (true) {
        idx = fastMod(idx, m_capacity);
        assert(idx < m_capacity);

        Cell const &cell = m_cells[idx];

        actual.key = cell.key;

        actual.value = cell.value;

        if (actual.key == key) {
            return actual.value;
        }

        if (actual.key == ParallelHashMap::ht_invalid_key) {
            return ParallelHashMap::ht_invalid_value;
        }

        idx++;
    }
}

bool HTAtomic128::update(uint64_t const key, uint64_t const value) {
    Cell const desired = makeCell(key, value);
    Cell expected = invalidCell();
    uint64_t const hashed_key = hash64(key);
    size_t idx = scaling_log(hashed_key, m_log_capacity);

    while (true) {
        idx = fastMod(idx, m_capacity);
        assert(idx < m_capacity);

        Cell &cell = m_cells[idx];

#pragma omp atomic read
        expected.key = cell.key;

#pragma omp atomic read
        expected.value = cell.value;

        if (expected.key == desired.key) {
            // If DWCAS fails we simply try again for the same cell.
            if (dwcas(cell, expected, desired)) {
                return true;
            }
        } else if (expected.key == ParallelHashMap::ht_invalid_key) {
            return false;
        } else {
            idx++;
        }
    }
}

HTAtomic128::Iterator::Iterator(Cells const &cells)
    : m_cells_it(std::begin(cells)), m_cells_end(std::end(cells)) {
    while (m_cells_it != m_cells_end && m_cells_it->key == ParallelHashMap::ht_invalid_key) {
        m_cells_it++;
    }
}

HTAtomic128::Iterator::pointer HTAtomic128::Iterator::operator->() {
    return &*m_cells_it;
}

HTAtomic128::Iterator &HTAtomic128::Iterator::operator++() {
    next();
    return *this;
}

HTAtomic128::Iterator HTAtomic128::Iterator::operator++(int) {
    Iterator tmp = *this;
    next();
    return tmp;
}

void HTAtomic128::Iterator::invalidate() {
    m_cells_it = m_cells_end;
}

void HTAtomic128::Iterator::next() {
    if (m_cells_it == m_cells_end) {
        return;
    }

    do {
        m_cells_it++;
    } while (m_cells_it != m_cells_end && m_cells_it->key == ParallelHashMap::ht_invalid_key);
}

HTAtomic128::Iterator HTAtomic128::end() const {
    HTAtomic128::Iterator end_it{m_cells};
    end_it.invalidate();
    return end_it;
}

size_t HTAtomic128::incrementGlobalOccupancy(size_t increment) {
    return Aux::incrementAtomically(m_global_occupancy, increment);
}

size_t HTAtomic128::globalOccupancy() const {
    return m_global_occupancy.load();
}

size_t HTAtomic128::capacity() const {
    return m_capacity;
}

bool HTAtomic128::filled() const {
    return ht_filled(m_global_occupancy.load(), m_capacity);
}

void HTAtomic128::roam(HTAtomic128 &target, uint32_t const p_count, uint32_t const p_id) {
    auto range_to_move = clusterRange(p_count, p_id);
    moveCells(range_to_move.first, range_to_move.second, target);
}

void HTAtomic128::moveCells(HTAtomic128::Cells::const_iterator begin,
                            HTAtomic128::Cells::const_iterator end, HTAtomic128 &target) {
    size_t local_occupancy = 0;
    for (HTAtomic128::Cells::const_iterator c = begin; c != end; ++c) {
        if (c->key != ParallelHashMap::ht_invalid_key) {
            bool inserted = target.insert(c->key, c->value);
            assert(inserted);
            local_occupancy += size_t(inserted);
        }
    }
    target.incrementGlobalOccupancy(local_occupancy);
}

std::pair<ParallelHashMap::HTAtomic128::Cells::const_iterator,
          ParallelHashMap::HTAtomic128::Cells::const_iterator>
HTAtomic128::clusterRange(uint32_t const p_count, uint32_t const p_id) {
    auto cells_begin = m_cells.begin();
    auto cells_end = m_cells.end();

    if (cells_begin == cells_end) {
        return std::make_pair(cells_end, cells_end);
    }

    HTAtomic128::Cells::difference_type const cell_count = std::distance(cells_begin, cells_end);
    if (cell_count < (p_count * p_count)) {
        if (p_id == 0) {
            return std::make_pair(cells_begin, cells_end);
        }

        return std::make_pair(cells_end, cells_end);
    }

    size_t const cells_per_thread = cell_count / p_count;

    auto p_cells_end = cells_begin;
    if (p_id == (p_count - 1)) {
        p_cells_end = cells_end;
    } else {
        size_t const max_advance = std::min<size_t>((p_id + 1) * cells_per_thread, cell_count);
        std::advance(p_cells_end, max_advance);

        while (p_cells_end != cells_end && p_cells_end->key != ParallelHashMap::ht_invalid_key) {
            ++p_cells_end;
        }
    }

    auto p_cells_begin = cells_begin;
    if (p_id > 0) {
        size_t const max_advance =
            std::min<size_t>(p_id * cells_per_thread, std::distance(p_cells_begin, p_cells_end));
        std::advance(p_cells_begin, max_advance);

        while (p_cells_begin != p_cells_end
               && p_cells_begin->key != ParallelHashMap::ht_invalid_key) {
            ++p_cells_begin;
        }
    }

    return std::make_pair(p_cells_begin, p_cells_end);
}

HTSyncData::HTSyncData(std::unique_ptr<HTAtomic128> &_source, std::unique_ptr<HTAtomic128> &_target,
                       std::atomic_uint32_t &_busy, std::atomic_bool &_request_growth, int _p_count,
                       int _p_id, uint32_t _insert_threshold)
    : source(_source), target(_target), busy(_busy), insert_threshold(_insert_threshold),
      request_growth(_request_growth), p_count(_p_count), p_id(_p_id) {}

HTHandle::HTHandle(HTSyncData sync_data) : m_ht(sync_data.source.get()), m_sync_data(sync_data) {
    setBitAtomically(m_sync_data.busy, m_sync_data.p_id);
}

HTHandle::~HTHandle() {
    size_t const occupancy = m_ht->incrementGlobalOccupancy(m_sync_data.insert_counter);

    if (ht_filled(occupancy, m_ht->capacity())) {
        m_ht = grow_hashtable(m_sync_data.source, m_sync_data.target, m_sync_data.request_growth,
                              m_sync_data.p_count, m_sync_data.p_id);
    }

    Aux::unsetBitAtomically(m_sync_data.busy, m_sync_data.p_id);

    while (m_sync_data.busy.load() != 0u) {
        if (m_sync_data.request_growth.load()) {
            grow_hashtable(m_sync_data.source, m_sync_data.target, m_sync_data.request_growth,
                           m_sync_data.p_count, m_sync_data.p_id);
        }
    }
}

bool HTHandle::insert(uint64_t const key, uint64_t const value) {
    bool const success = m_ht->insert(key, value);
    m_sync_data.insert_counter += uint32_t(success);

    if (m_sync_data.insert_threshold == m_sync_data.insert_counter) {
        size_t const occupancy = m_ht->incrementGlobalOccupancy(m_sync_data.insert_counter);

        if (ht_filled(occupancy, m_ht->capacity())) {
            m_ht =
                grow_hashtable(m_sync_data.source, m_sync_data.target, m_sync_data.request_growth,
                               m_sync_data.p_count, m_sync_data.p_id);
        }

        m_sync_data.insert_counter = 0;
    }

    return success;
}

uint64_t HTHandle::find(uint64_t const key) const {
    return m_ht->find(key);
}

bool HTHandle::update(uint64_t const key, uint64_t const value) const {
    return m_ht->update(key, value);
}

HTAtomic128 const &HTHandle::hashtable() {
    return *m_ht;
}

ParallelHashMap::ParallelHashMap() : ParallelHashMap(ParallelHashMap::ht_begin_capacity) {}

ParallelHashMap::ParallelHashMap(size_t const begin_capacity) {
    if (begin_capacity == 0 || ((begin_capacity & (begin_capacity - 1)) != 0)) {
        throw std::runtime_error("begin_capacity must be a power of 2 and greater than 0");
    }
    m_source = std::make_unique<HTAtomic128>(begin_capacity);
}

ParallelHashMap::ParallelHashMap(ParallelHashMap const &other) {
    m_source = std::make_unique<HTAtomic128>(*other.m_source.get());
}

ParallelHashMap::ParallelHashMap(ParallelHashMap &&other) noexcept : ParallelHashMap() {
    swap(*this, other);
}

HTHandle ParallelHashMap::makeHandle() {
    Aux::HTSyncData sync_data{m_source,
                              m_target,
                              m_busy_bitset,
                              m_request_growth,
                              omp_get_num_threads(),
                              omp_get_thread_num(),
                              randomThreadRange()};

    // auto handle = std::make_unique<HTHandle>(sync_data);
    HTHandle handle(std::move(sync_data));
    return handle;
    // return handle;
}

HTAtomic128 const *ParallelHashMap::currentTable() const {
    return m_source.get();
}

uint32_t ParallelHashMap::randomThreadRange() {
    std::uniform_int_distribution<uint32_t> m_dist(1, omp_get_num_threads());
    return m_dist(m_generator);
}

} // namespace Aux
