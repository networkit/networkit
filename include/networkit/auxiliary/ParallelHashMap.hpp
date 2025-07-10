#ifndef NETWORKIT_AUXILIARY_PARALLEL_HASH_MAP_HPP_
#define NETWORKIT_AUXILIARY_PARALLEL_HASH_MAP_HPP_

#include <omp.h>

#include <limits>
#include <memory>
#include <random>
#include <vector>

#include <networkit/auxiliary/AtomicUtils.hpp>
#include <networkit/auxiliary/HashUtils.hpp>

#if defined(__GNUC__)
#define ASM_CLOBBER_CC
#else
#define ASM_CLOBBER_CC "cc"
#endif

namespace Aux {

// An object of this class will handle a dynamic hashtable and offers factory
// methods to acquire hashtable resources.
class ParallelHashMap {
public:
    friend void swap(ParallelHashMap &p, ParallelHashMap &q) noexcept {
        using std::swap;

        swap(p.m_source, q.m_source);
        swap(p.m_target, q.m_target);
        swapAtomicsNonAtomically<uint32_t>(p.m_busy_bitset, q.m_busy_bitset);
        swapAtomicsNonAtomically<bool>(p.m_request_growth, q.m_request_growth);
    }

    class HTAtomic128;
    class HTSyncData;
    class HTHandle;

    static constexpr uint64_t ht_invalid_key = std::numeric_limits<uint64_t>::max();
    static constexpr uint64_t ht_invalid_value = std::numeric_limits<uint64_t>::max() - 1;

    // This value shall represent some storage capacity that can be handled within a
    // CPU cache (e.g. private MLC).
    static constexpr size_t const ht_begin_capacity = 4096u;

    static constexpr uint64_t ht_key_space = sizeof(uint64_t) * 8;

    /** Default constructor */
    ParallelHashMap();

    /**
     *  Constructs a ParallelHashMap with a given initial capacity.
     *  The capacity must be a power of 2 and greater than 0.
     *  @param begin_capacity The initial capacity of the hashtable.
     *  @throws std::runtime_error if begin_capacity is not a power of 2 or
     *          if begin_capacity is 0.
     * @note The initial capacity should also exceed the square number of threads. For example,
     *       if you have 4 threads, the initial capacity should be at least 32 (4^2 = 16, next power
     * of 2 is 32). If the initial capacity is too small, the hashtable might get locked during
     * inserts. The default value of 4096 (64^2) is a good starting point for most applications.
     */
    ParallelHashMap(size_t begin_capacity);

    ~ParallelHashMap() = default;

    ParallelHashMap(ParallelHashMap const &other);

    ParallelHashMap(ParallelHashMap &&other) noexcept;

    ParallelHashMap &operator=(ParallelHashMap other) {
        swap(*this, other);
        return *this;
    }

    HTHandle makeHandle();
    /**
     *  Returns the current hashtable that is being used by the handle.
     *  This is the source table that is being used for inserts and lookups.
     *  @return A pointer to the current hashtable.
     *  @note This pointer is not guaranteed to be valid after the handle has been
     *        created, as the hashtable may be swapped out during operations.
     */
    HTAtomic128 const *currentTable() const;

private:
    std::unique_ptr<HTAtomic128> m_source;
    std::unique_ptr<HTAtomic128> m_target;
    std::atomic_uint32_t m_busy_bitset{0u};
    std::atomic_bool m_request_growth{false};

    // Used to determine how often to check for the hashtables fill factor for
    // every 1 to p inserts. We're keeping this an internal random number
    // generator engine. However, this could also be injected when constructing.
    std::mt19937 m_generator;

    uint32_t randomThreadRange();
};

// This is the Hashtable Atomic 128 (HTAtomic128) class which stores (key,
// value) pairs, both of 64 bit and operates on them atomically using DWCAS
// operations.
//
// The basis of this implementation is the "Concurrent Hash Tables: Fast and
// general(?)!" from Maier et al. (2016).
//
// Some interfaces are not needed for our case such as insertOrUpdate() or
// update(). We expect the key's never to be updated during runtime: the handle
// for vertex u will never change during the execution of the program.
class ParallelHashMap::HTAtomic128 {
public:
    // In order to use DWCAS we must align our Cell within a 16 bytes window.
    struct alignas(16) Cell {
        uint64_t key{ParallelHashMap::ht_invalid_key};
        uint64_t value{ParallelHashMap::ht_invalid_value};
    };

    using Cells = std::vector<Cell>;

    HTAtomic128();
    HTAtomic128(size_t capacity);
    ~HTAtomic128() = default;
    HTAtomic128(HTAtomic128 const &other);
    HTAtomic128(HTAtomic128 &&other) noexcept;

    HTAtomic128 &operator=(HTAtomic128 other);

    // Thanks to GManNickG on Stackoverflow explaining the copy-swap idiom:
    // https://stackoverflow.com/a/3279550
    friend void swap(HTAtomic128 &a, HTAtomic128 &b) noexcept {
        using std::swap;

        swap(a.m_cells, b.m_cells);
        swap(a.m_capacity, b.m_capacity);
        swap(a.m_log_capacity, b.m_log_capacity);
        swapAtomicsNonAtomically(a.m_global_occupancy, b.m_global_occupancy);
    }

    /** Inserts a key-value pair into the hashtable.
     * If a new entry could not be created, the insertion will fail and return false.
     * This means that the key already exists in the hashtable. Use update() to change the value
     * associated with the key.
     * @param key The key to insert.
     * @param value The value to insert.
     * @return True if the insertion was successful, false if the key already exists.
     * @note This function is thread-safe and can be called concurrently.
     */
    bool insert(uint64_t key, uint64_t value);

    /**
     * Finds the value associated with a key in the hashtable.
     * If the key is not present, it returns ht_invalid_key.
     * @param key The key to find.
     * @return The value associated with the key, or ht_invalid_key if the key is not found.
     * @note This function is thread-safe and can be called concurrently.
     */
    uint64_t find(uint64_t key) const;

    /**
     * Updates the value associated with a key in the hashtable.
     * If the key does not exist, the update will fail and return false.
     * If the key exists, the value will be updated to the new value.
     * @param key The key to update.
     * @param value The new value to associate with the key.
     * @return True if the update was successful, false if the key does not exist.
     * @note This function is thread-safe and can be called concurrently.
     */
    bool update(uint64_t key, uint64_t value);

    struct Iterator {
        using iterator_category = std::input_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Cell const;
        using pointer = Cell const *;
        using reference = Cell const &;

        Iterator(Cells const &cells);

        friend void swap(Iterator &a, Iterator &b) {
            using std::swap;

            swap(a.m_cells_it, b.m_cells_it);
            swap(a.m_cells_end, b.m_cells_end);
        }

        reference operator*() { return *m_cells_it; }
        pointer operator->();

        // Prefix increment operator
        Iterator &operator++();

        // Postfix increment operator
        Iterator operator++(int);

        bool operator==(const HTAtomic128::Iterator &b) const {
            return this->m_cells_it == b.m_cells_it;
        }

        bool operator!=(const HTAtomic128::Iterator &b) const { return !(*this == b); }

        void invalidate();

    private:
        void next();

    private:
        Cells::const_iterator m_cells_it;
        Cells::const_iterator m_cells_end;
    };

    Iterator begin() const { return HTAtomic128::Iterator{m_cells}; }

    Iterator end() const;

    /**
     *  Increments the global occupancy counter by the given amount.
     *  @param increment The amount to increment the counter by.
     *  @return The new value of the global occupancy counter.
     */
    size_t incrementGlobalOccupancy(size_t increment);
    /**
     *  Checks if the hashtable is filled based on the occupancy and capacity.
     *  @param occupancy The current occupancy of the hashtable.
     *  @param capacity The capacity of the hashtable.
     *  @return true if the hashtable is filled, false otherwise.
     */
    size_t globalOccupancy() const;

    /**
     * Returns the current capacity of the hashtable.
     * The capacity is the number of cells in the hashtable.
     * It is always a power of 2.
     * @return The current capacity of the hashtable.
     */
    size_t capacity() const;

    /**
     *  Returns the cells of the hashtable.
     *  This is a reference to the internal vector of cells.
     *  @return A reference to the vector of cells.
     */
    Cells const &cells() const { return m_cells; }

    /**
     * Returns true if the hashtable is filled, meaning that the occupancy
     * reaches greater equal than half of the capacity (alpha = 0.5 * capacity).
     * @return True if the hashtable is filled, false otherwise.
     */
    bool filled() const;

    /**
     *  Roams the source hashtable to the target hashtable.
     *  This operation is performed in parallel, with each thread
     *  operating on a defined space of the source table.
     *  The target table must be empty before calling this function.
     *  @param target The target hashtable to roam to.
     *  @param p_count The number of threads in the parallel region.
     *  @param p_id The id of the calling thread.
     */
    void roam(ParallelHashMap::HTAtomic128 &target, uint32_t p_count = 1, uint32_t p_id = 0);

    /**
     *  Returns a range of cells that belong to the cluster of the given thread id.
     *  The range is defined by the number of threads and the thread id.
     *  The range is guaranteed to be valid and non-empty for the first thread (p_id = 0).
     *  For other threads, the range may be empty if there are not enough cells to
     *  fill the cluster.
     *  @param p_count The number of threads in the parallel region.
     *  @param p_id The id of the thread for which the cluster range is requested.
     *  @return A pair of iterators defining the range of cells for the cluster.
     */
    std::pair<ParallelHashMap::HTAtomic128::Cells::const_iterator,
              ParallelHashMap::HTAtomic128::Cells::const_iterator>
    clusterRange(uint32_t p_count = 1, uint32_t p_id = 0);

private:
    static Cell makeCell(uint64_t key, uint64_t value) {
        Cell c;
        c.key = key;
        c.value = value;
        return c;
    }

    static Cell invalidCell() {
        return makeCell(ParallelHashMap::ht_invalid_key, ParallelHashMap::ht_invalid_value);
    }

    // Double Width Compare And Swap:
    //
    // If r matches expected, write desired to r and return true. Otherwise,
    // read r into expected and return false.
    //
    // We using this user defined dwcas since we can't be sure we're operating
    // atomically on two 64-bit values: https://timur.audio/dwcas-in-c
    static bool dwcas(HTAtomic128::Cell &r, HTAtomic128::Cell &expected,
                      HTAtomic128::Cell const &desired) {
        bool res;
#if __x86_64__
        asm volatile("lock cmpxchg16b %1"
                     : "=@ccz"(res), "+m"(r), "+a"(expected.key), "+d"(expected.value)
                     : "b"(desired.key), "c"(desired.value)
                     : "memory");
#elif defined(_MSC_VER)
        res = _InterlockedCompareExchange128(
            reinterpret_cast<volatile __int64 *>(&r), static_cast<__int64>(desired.value),
            static_cast<__int64>(desired.key), reinterpret_cast<__int64 *>(&expected));

#else
        // For ARM64 we use the ld/stxp instruction pair to perform the DWCAS operation.
        // See:
        // https://github.com/boostorg/atomic/blob/master/include/boost/atomic/detail/core_arch_ops_gcc_aarch64.hpp
        HTAtomic128::Cell original;
        asm volatile("mov %w[success], #0\n\t"
                     "ldxp %x[original_0], %x[original_1], %[storage]\n\t"
                     "cmp %x[original_0], %x[expected_0]\n\t"
                     "ccmp %x[original_1], %x[expected_1], #0, eq\n\t"
                     "b.ne 1f\n\t"
                     "stxp %w[success], %x[desired_0], %x[desired_1], %[storage]\n\t"
                     "eor %w[success], %w[success], #1\n\t"
                     "1:\n\t"
                     : [success] "=&r"(res), [storage] "+Q"(r), [original_0] "=&r"(original.key),
                       [original_1] "=&r"(original.value)
                     : [desired_0] "r"(desired.key), [desired_1] "r"(desired.value),
                       [expected_0] "r"(expected.key), [expected_1] "r"(expected.value)
                     : ASM_CLOBBER_CC "memory");
#endif
        return res;
    }

    void moveCells(ParallelHashMap::HTAtomic128::Cells::const_iterator begin,
                   ParallelHashMap::HTAtomic128::Cells::const_iterator end,
                   ParallelHashMap::HTAtomic128 &target);

private:
    size_t m_capacity{ParallelHashMap::ht_begin_capacity};
    uint64_t m_log_capacity{static_cast<uint64_t>(std::log2(m_capacity))};
    Cells m_cells;
    std::atomic_size_t m_global_occupancy{0};
};

class ParallelHashMap::HTSyncData {
public:
    HTSyncData(std::unique_ptr<HTAtomic128> &_source, std::unique_ptr<HTAtomic128> &_target,
               std::atomic_uint32_t &_busy, std::atomic_bool &_request_growth, int _p_count,
               int _p_id, uint32_t _insert_threshold);

    std::unique_ptr<HTAtomic128> &source;
    std::unique_ptr<HTAtomic128> &target;

    // Represents a bitset where the bit is set to 1 at position p_id if the
    // thread is still busy operating on hash table operations.
    std::atomic_uint32_t &busy;
    uint32_t insert_counter{0};
    uint32_t const insert_threshold;

    // Whenever this is set, the hashtable has to be grown by all threads of the
    // omp team.
    std::atomic_bool &request_growth;
    int const p_count;
    int const p_id;
};

class ParallelHashMap::HTHandle {
public:
    HTHandle(HTSyncData sync_data);

    ~HTHandle();

    bool insert(uint64_t key, uint64_t value);

    uint64_t find(uint64_t key) const;

    bool update(uint64_t key, uint64_t value) const;

    HTAtomic128 const &hashtable();

private:
    HTAtomic128 *m_ht;
    HTSyncData m_sync_data;
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_PARALLEL_HASH_MAP_HPP_
