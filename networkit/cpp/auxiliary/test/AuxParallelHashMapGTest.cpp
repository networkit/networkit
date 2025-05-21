/*
 * AuxGTest.cpp
 *
 *  Created on: 20.05.2025
 *      Author: Fabian Brandt-Tumescheit
 *              Florian Willich
 */

#include <gtest/gtest.h>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/AtomicUtils.hpp>
#include <networkit/auxiliary/ParallelHashMap.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

using MockupData = std::vector<std::pair<node, size_t>>;
using HTAtomic128 = Aux::ParallelHashMap::HTAtomic128;
using HTHandle = Aux::ParallelHashMap::HTHandle;
using HTSyncData = Aux::ParallelHashMap::HTSyncData;

class AuxParallelGrowingHTGTest : public testing::Test {

public:
    bool get_and_check_all(std::vector<node> const &mockup_data, HTAtomic128 &ht) {
        size_t idx = 0;
        bool all_good = true;
        for (size_t i = 0; i < mockup_data.size() && all_good; ++i) {
            idx = ht.find(mockup_data[i]);

            all_good = idx < mockup_data.size();
            all_good = all_good && mockup_data[idx] == mockup_data[i];
        }

        return all_good;
    }

    MockupData generateMockupData(size_t const fill_64bit_values = 148) {
        MockupData mockup_data;

        constexpr node minimum_of_range = 1;
        constexpr node maximum_of_range = 356;
        auto &prng = Aux::Random::getURNG();
        std::uniform_int_distribution<uint64_t> distr{minimum_of_range, maximum_of_range};

        node u = 0;
        for (size_t i = 0; i < fill_64bit_values; ++i) {
            u += distr(prng);
            mockup_data.push_back(std::make_pair(u, i));
        }

        return mockup_data;
    }

    bool checkEntriesForMockupData(MockupData const &mockup_data, HTAtomic128 const &ht) {
        bool all_good = true;
        for (size_t i = 0; i < mockup_data.size() && all_good; ++i) {
            size_t const retrieved_idx = ht.find(mockup_data[i].first);

            all_good = retrieved_idx < mockup_data.size();
            all_good = all_good && retrieved_idx == mockup_data[i].second;
        }
        return all_good;
    }

    bool insert_data(HTAtomic128 &ht, MockupData const &data) {
        bool all_inserted_successful = true;
        for (size_t i = 0; i < data.size(); ++i) {
            all_inserted_successful = ht.insert(data[i].first, data[i].second);
            if (!all_inserted_successful) {
                break;
            }
        }
        return all_inserted_successful;
    }

    bool insert_data(HTHandle &handle, MockupData const &data) {
        bool all_inserted_successful = true;
        for (size_t i = 0; i < data.size() && all_inserted_successful; ++i) {
            all_inserted_successful = handle.insert(data[i].first, data[i].second);
        }
        return all_inserted_successful;
    }
};

TEST_F(AuxParallelGrowingHTGTest, testHashUtilsfastMod) {
    uint32_t idx = 8;
    size_t capacity = 16;

    uint32_t idx_wrapped = Aux::fast_mod(idx, capacity);

    ASSERT_EQ(idx_wrapped, idx);

    idx = 13;
    idx_wrapped = Aux::fast_mod(idx, capacity);

    ASSERT_EQ(idx_wrapped, idx);

    idx = 17;
    idx_wrapped = Aux::fast_mod(idx, capacity);

    ASSERT_EQ(idx_wrapped, 1);
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128Capacity) {
    MockupData const mockup_data = generateMockupData(125000);

    HTAtomic128 ht{131072};

    for (size_t i = 0; i < mockup_data.size(); ++i) {
        ht.insert(mockup_data[i].first, mockup_data[i].second);
    }

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128Empty) {
    HTAtomic128 ht{Aux::ParallelHashMap::ht_begin_capacity};

    ASSERT_EQ(ht.find(0), Aux::ParallelHashMap::ht_invalid_value);
    ASSERT_EQ(ht.find(155), Aux::ParallelHashMap::ht_invalid_value);
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128Constructor) {
    std::vector<node> mockupData;
    constexpr size_t fill_value_with_64bit_values = 148;
    constexpr node vertex_offset_from_index = 25;

    for (size_t i = 0; i < fill_value_with_64bit_values; ++i) {
        mockupData.push_back(i + vertex_offset_from_index);
    }

    HTAtomic128 ht{Aux::ParallelHashMap::ht_begin_capacity};

    for (size_t i = 0; i < mockupData.size(); ++i) {
        ht.insert(mockupData[i], i);
    }

    ASSERT_TRUE(get_and_check_all(mockupData, ht));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128CopyConstructor) {
    std::vector<node> mockupData;
    constexpr size_t fill_value_with_64bit_values = 148;
    constexpr node vertex_offset_from_index = 25;

    for (size_t i = 0; i < fill_value_with_64bit_values; ++i) {
        mockupData.push_back(i + vertex_offset_from_index);
    }

    HTAtomic128 ht{Aux::ParallelHashMap::ht_begin_capacity};

    for (size_t i = 0; i < mockupData.size(); ++i) {
        ht.insert(mockupData[i], i);
    }

    HTAtomic128 dvh_copy(ht);

    ASSERT_TRUE(get_and_check_all(mockupData, dvh_copy));
    ASSERT_TRUE(get_and_check_all(mockupData, ht));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128MoveConstructor) {
    std::vector<node> mockupData;
    constexpr size_t fill_value_with_64bit_values = 148;
    constexpr node vertex_offset_from_index = 25;

    for (size_t i = 0; i < fill_value_with_64bit_values; ++i) {
        mockupData.push_back(i + vertex_offset_from_index);
    }

    HTAtomic128 ht{Aux::ParallelHashMap::ht_begin_capacity};

    for (size_t i = 0; i < mockupData.size(); ++i) {
        ht.insert(mockupData[i], i);
    }

    HTAtomic128 dvh_moved(std::move(ht));

    ASSERT_TRUE(get_and_check_all(mockupData, dvh_moved));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128Swap) {
    std::vector<node> mockupData;
    constexpr size_t fill_value_with_64bit_values = 148;
    constexpr node vertex_offset_from_index = 25;

    for (size_t i = 0; i < fill_value_with_64bit_values; ++i) {
        mockupData.push_back(i + vertex_offset_from_index);
    }

    HTAtomic128 ht{Aux::ParallelHashMap::ht_begin_capacity};

    for (size_t i = 0; i < mockupData.size(); ++i) {
        ht.insert(mockupData[i], i);
    }

    std::vector<node> mockupDataB;
    constexpr size_t fill_value_with_64bit_values_b = 54;
    constexpr node vertex_offset_from_index_b = 26;

    for (size_t i = 0; i < fill_value_with_64bit_values_b; ++i) {
        mockupDataB.push_back(i + vertex_offset_from_index_b);
    }

    HTAtomic128 ht_b{Aux::ParallelHashMap::ht_begin_capacity};

    for (size_t i = 0; i < mockupDataB.size(); ++i) {
        ht_b.insert(mockupDataB[i], i);
    }

    ht_b.increment_global_occupancy(mockupDataB.size());

    swap(ht, ht_b);

    ASSERT_TRUE(get_and_check_all(mockupDataB, ht));
    ASSERT_TRUE(get_and_check_all(mockupData, ht_b));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128CopyAssignment) {
    std::vector<node> mockupData;
    constexpr size_t fill_value_with_64bit_values = 148;
    constexpr node vertex_offset_from_index = 25;

    for (size_t i = 0; i < fill_value_with_64bit_values; ++i) {
        mockupData.push_back(i + vertex_offset_from_index);
    }

    HTAtomic128 ht{Aux::ParallelHashMap::ht_begin_capacity};

    for (size_t i = 0; i < mockupData.size(); ++i) {
        ht.insert(mockupData[i], i);
    }
    HTAtomic128 ht_b{Aux::ParallelHashMap::ht_begin_capacity};
    ht_b = ht;

    ASSERT_TRUE(get_and_check_all(mockupData, ht_b));
    ASSERT_TRUE(get_and_check_all(mockupData, ht));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128InputIterator) {
    MockupData const mockup_data = generateMockupData();

    HTAtomic128 ht{Aux::ParallelHashMap::ht_begin_capacity};
    bool all_inserted_successful = true;
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = ht.insert(mockup_data[i].first, mockup_data[i].second);
    }

    ASSERT_TRUE(all_inserted_successful);

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht));

    // Begin and End Iterator
    HTAtomic128::Iterator it = ht.begin();
    HTAtomic128::Iterator it_end = ht.end();

    ASSERT_TRUE(it != it_end);

    auto it_2 = ht.begin();

    ASSERT_TRUE(it == it_2);

    auto it_end_2 = ht.end();

    ASSERT_TRUE(it_end == it_end_2);

    it.invalidate();

    ASSERT_TRUE(it == it_end);

    // Manual check
    it = ht.begin();
    size_t iterated_elements = 0;

    bool all_valid = true;
    bool not_end = true;
    while (iterated_elements < mockup_data.size() && all_valid && not_end) {
        all_valid = it->key != Aux::ParallelHashMap::ht_invalid_key;
        not_end = it != ht.end();
        it++;
        iterated_elements++;
    }

    ASSERT_TRUE(all_valid);
    ASSERT_TRUE(not_end);
    ASSERT_EQ(iterated_elements, mockup_data.size());

    it++;
    ASSERT_TRUE(it == ht.end());

    // Iterate
    bool keys_are_correct = true;
    bool values_are_correct = true;
    size_t elements_iterated = 0;

    bool it_not_end = true;
    for (auto it = ht.begin(); it != ht.end() && keys_are_correct && values_are_correct; ++it) {
        keys_are_correct = it->key == mockup_data[it->value].first;
        values_are_correct = it->value == mockup_data[it->value].second;
        it_not_end = it != ht.end();
        elements_iterated++;
    }

    ASSERT_TRUE(it_not_end);
    ASSERT_TRUE(keys_are_correct);
    ASSERT_TRUE(values_are_correct);
    ASSERT_EQ(elements_iterated, mockup_data.size());

    // For
    elements_iterated = 0;

    keys_are_correct = true;
    values_are_correct = true;
    for (auto c : ht) {
        keys_are_correct = keys_are_correct && (c.key == mockup_data[c.value].first);
        values_are_correct = values_are_correct && (c.value == mockup_data[c.value].second);
        elements_iterated++;
    }

    ASSERT_TRUE(keys_are_correct);
    ASSERT_TRUE(values_are_correct);
    ASSERT_TRUE(elements_iterated == mockup_data.size());
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128ConcurrentAccess) {
    Aux::setNumberOfThreads(2);

    constexpr size_t one_megabyte_filled_with_64bit_values = 125000;

    auto mockup_data = generateMockupData(one_megabyte_filled_with_64bit_values);

    constexpr size_t begin_capacity = 131072;
    HTAtomic128 ht{begin_capacity};

#pragma omp parallel for
    for (auto e : mockup_data) {
        ht.insert(e.first, e.second);
    }

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128ConcurrentAccessWithRandomWrites) {
    Aux::setNumberOfThreads(2);

    HTAtomic128 ht_r;

#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        if (thread_id == 0) {
            ht_r.insert(0, 0);
            ht_r.insert(1, 0);
        }

        if (thread_id == 1) {
            ht_r.insert(0, 1);
            ht_r.insert(1, 1);
        }
    }

    node zero = ht_r.find(0);
    ASSERT_TRUE((zero == 0 || zero == 1));
    node one = ht_r.find(1);
    ASSERT_TRUE((one == 0 || one == 1));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128ConcurrentAccessWithBarrier) {
    Aux::setNumberOfThreads(2);

    HTAtomic128 ht_l;
    uint64_t foundZero = 0;
    uint64_t foundOne = 0;

#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        if (thread_id == 0) {
            ht_l.insert(0, 0);
            ht_l.insert(1, 0);
        }

#pragma omp barrier
        if (thread_id == 1) {
            foundZero = ht_l.find(0);
            foundOne = ht_l.find(1);
        }
    }

    ASSERT_EQ(foundZero, 0);
    ASSERT_EQ(foundOne, 0);
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128ConcurrentAccessWithUpdate) {
    Aux::setNumberOfThreads(2);

    HTAtomic128 ht_u;

    ht_u.insert(0, 0);
    ht_u.insert(1, 0);

#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        if (thread_id == 0) {
            ht_u.update(0, 2);
        } else {
            ht_u.update(0, 3);
        }
    }

    uint64_t const updated_value = ht_u.find(0);

    ASSERT_TRUE((updated_value == 2 || updated_value == 3));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128InsertDualThread) {
    Aux::setNumberOfThreads(2);

    constexpr size_t fill_64bit_values = 10254;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

    constexpr size_t begin_capacity = 16384u;
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = true;
#pragma omp parallel
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = ht_first.insert(mockup_data[i].first, mockup_data[i].second);
    }

    ASSERT_TRUE(all_inserted_successful);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_first));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128ClusterRangeSingleThread) {
    constexpr size_t fill_64bit_values = 345;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

    constexpr size_t begin_capacity = 512;
    HTAtomic128 source{begin_capacity};

    bool all_inserted_successful = true;
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = source.insert(mockup_data[i].first, mockup_data[i].second);
    }
    ASSERT_TRUE(all_inserted_successful);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, source));

    auto const cells_begin = std::cbegin(source.cells());
    auto const cells_end = std::cend(source.cells());

    auto c_range = source.cluster_range();
    ASSERT_EQ(c_range.first, cells_begin);
    ASSERT_EQ(c_range.second, cells_end);
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128ClusterRangeDualThread) {
    constexpr size_t fill_64bit_values = 345;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

    constexpr size_t begin_capacity = 512;
    HTAtomic128 source{begin_capacity};

    bool all_inserted_successful = true;
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = source.insert(mockup_data[i].first, mockup_data[i].second);
    }
    ASSERT_TRUE(all_inserted_successful);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, source));

    auto const cells_begin = std::cbegin(source.cells());
    auto const cells_end = std::cend(source.cells());

    auto c_range_t0 = source.cluster_range(2, 0);
    ASSERT_EQ(c_range_t0.first, cells_begin);
    ASSERT_TRUE(c_range_t0.second != cells_end);

    auto c_range_t1 = source.cluster_range(2, 1);
    ASSERT_TRUE(c_range_t1.first != cells_begin);
    ASSERT_EQ(c_range_t1.second, cells_end);

    ASSERT_EQ(c_range_t0.second, c_range_t1.first);
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128ClusterRangeQuadrupleThread) {
    constexpr size_t fill_64bit_values = 345;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

    constexpr size_t begin_capacity = 512;
    HTAtomic128 source{begin_capacity};

    bool all_inserted_successful = true;
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = source.insert(mockup_data[i].first, mockup_data[i].second);
    }
    ASSERT_TRUE(all_inserted_successful);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, source));

    auto const cells_begin = std::cbegin(source.cells());
    auto const cells_end = std::cend(source.cells());

    uint32_t const thread_count = 4;

    auto c_range_t0 = source.cluster_range(thread_count, 0);
    ASSERT_EQ(c_range_t0.first, cells_begin);
    ASSERT_TRUE(c_range_t0.second != cells_end);

    auto c_range_t1 = source.cluster_range(thread_count, 1);
    ASSERT_TRUE(c_range_t1.first != cells_begin);
    ASSERT_TRUE(c_range_t1.second != cells_end);
    ASSERT_EQ(c_range_t0.second, c_range_t1.first);

    auto c_range_t2 = source.cluster_range(thread_count, 2);
    ASSERT_TRUE(c_range_t2.first != cells_begin);
    ASSERT_TRUE(c_range_t2.second != cells_end);
    ASSERT_EQ(c_range_t1.second, c_range_t2.first);

    auto c_range_t3 = source.cluster_range(thread_count, 3);
    ASSERT_TRUE(c_range_t3.first != cells_begin);
    ASSERT_EQ(c_range_t3.second, cells_end);
    ASSERT_EQ(c_range_t2.second, c_range_t3.first);
}

GTEST_TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128RoamingSequential) {
    constexpr size_t fill_64bit_values = 23;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

    constexpr size_t begin_capacity = 32;
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = true;
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = ht_first.insert(mockup_data[i].first, mockup_data[i].second);
    }

    ASSERT_TRUE(all_inserted_successful);

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_first));

    constexpr size_t next_capacity = begin_capacity << 1;
    HTAtomic128 ht_second{next_capacity};

    ht_first.roam(ht_second);

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_second));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128RoamingDualThread) {
    Aux::setNumberOfThreads(2);

    constexpr size_t fill_64bit_values = 23;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

    constexpr size_t begin_capacity = 32;
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = insert_data(ht_first, mockup_data);

    ASSERT_TRUE(all_inserted_successful);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_first));

    constexpr size_t next_capacity = begin_capacity << 1;
    HTAtomic128 ht_second{next_capacity};

#pragma omp parallel
    { ht_first.roam(ht_second, 2, omp_get_thread_num()); }

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_second));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128BigRoaming) {
    Aux::setNumberOfThreads(4);

    constexpr size_t fill_64bit_values = 10254;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

    constexpr size_t begin_capacity = 16384u; //
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = insert_data(ht_first, mockup_data);

    ASSERT_TRUE(all_inserted_successful);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_first));

    constexpr size_t next_capacity = begin_capacity << 1;
    HTAtomic128 ht_second{next_capacity};

#pragma omp parallel
    { ht_first.roam(ht_second, 4, omp_get_thread_num()); }

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_second));
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128HTHandleHTRoamingSingleThread) {
    MockupData const mockup_data = generateMockupData();
    std::mt19937 generator;

    Aux::setNumberOfThreads(1);

    std::unique_ptr<HTAtomic128> source =
        std::make_unique<HTAtomic128>(Aux::ParallelHashMap::ht_begin_capacity);
    std::unique_ptr<HTAtomic128> target;

    std::atomic_uint32_t busy_bitset{0u};
    std::atomic_bool request_growth{false};

    HTSyncData sync_data{
        source, target, busy_bitset, request_growth, 1, Aux::getCurrentNumberOfThreads(), 2};
    HTHandle handle{sync_data};

    bool const all_inserted_successful = insert_data(handle, mockup_data);

    ASSERT_TRUE(all_inserted_successful);

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, handle.hashtable()));

    ASSERT_EQ(handle.hashtable().global_occupancy(), mockup_data.size());
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128HTHandleHTRoamingDualThread) {
    MockupData const mockup_data = generateMockupData();
    std::mt19937 generator;

    Aux::setNumberOfThreads(2);

    std::unique_ptr<HTAtomic128> source =
        std::make_unique<HTAtomic128>(Aux::ParallelHashMap::ht_begin_capacity);
    std::unique_ptr<HTAtomic128> target;

    std::atomic_uint32_t busy_bitset{0u};
    std::atomic_bool request_growth{false};
    std::vector<bool> successful_inserts({true, true});
#pragma omp parallel
    {
        int const thread_id = omp_get_thread_num();
        HTSyncData sync_data{
            source,    target, busy_bitset, request_growth, Aux::getCurrentNumberOfThreads(),
            thread_id, 2};
        HTHandle handle{sync_data};

        auto mockup_data_begin = std::begin(mockup_data);
        std::advance(mockup_data_begin, thread_id == 0 ? 0 : mockup_data.size() / 2);

        auto mockup_data_end = std::begin(mockup_data);
        std::advance(mockup_data_end, thread_id == 0 ? mockup_data.size() / 2 : mockup_data.size());

        for (auto it = mockup_data_begin; it != mockup_data_end; ++it) {
            successful_inserts[thread_id] =
                successful_inserts[thread_id] && handle.insert(it->first, it->second);
        }
    }
    ASSERT_TRUE(successful_inserts[0]);
    ASSERT_TRUE(successful_inserts[1]);

    HTAtomic128 &ht = *source.get();

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht));

    ASSERT_EQ(ht.global_occupancy(), mockup_data.size());
}

TEST_F(AuxParallelGrowingHTGTest, testHTAtomic128HTHandleHTRoamingDualThreadRoaming) {
    MockupData const mockup_data = generateMockupData();
    std::mt19937 generator;

    Aux::setNumberOfThreads(2);

    std::unique_ptr<HTAtomic128> source = std::make_unique<HTAtomic128>(32u);
    std::unique_ptr<HTAtomic128> target;

    size_t const mockup_data_share = mockup_data.size() / 2;

    std::atomic_uint32_t busy_bitset{0u};
    std::atomic_bool request_growth{false};
    std::vector<bool> successful_inserts({true, true});

#pragma omp parallel
    {
        int const thread_id = omp_get_thread_num();
        HTSyncData sync_data{source, target, busy_bitset, request_growth, 2, thread_id, 2};
        HTHandle handle{sync_data};

        auto begin = std::begin(mockup_data);
        std::advance(begin, mockup_data_share * thread_id);

        auto end = begin;

        if (thread_id == 1) {
            end = std::end(mockup_data);
        } else {
            std::advance(end, mockup_data_share);
        }

        bool all_inserted_successful = true;
        for (auto it = begin; it != end; it++) {
            successful_inserts[thread_id] =
                successful_inserts[thread_id] && handle.insert(it->first, it->second);
        }
    }
    ASSERT_TRUE(successful_inserts[0]);
    ASSERT_TRUE(successful_inserts[1]);

    HTAtomic128 &ht_2 = *source.get();

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, ht_2));
}

TEST_F(AuxParallelGrowingHTGTest, testHTCustodianConstructor) {
    Aux::setNumberOfThreads(1);

    // Swap

    Aux::ParallelHashMap phm_a{32};
    Aux::ParallelHashMap phm_b{32};

    MockupData mockup_data = generateMockupData();

    auto handle_a = phm_a.make_handle();
    auto handle_b = phm_b.make_handle();

    ASSERT_TRUE(handle_a->insert(mockup_data[0].first, mockup_data[0].second));
    ASSERT_TRUE(handle_b->insert(mockup_data[1].first, mockup_data[1].second));

    swap(phm_a, phm_b);

    auto handle_a_new = phm_a.make_handle();
    auto handle_b_new = phm_b.make_handle();

    ASSERT_TRUE(handle_a_new->find(mockup_data[1].first) == mockup_data[1].second);
    ASSERT_TRUE(handle_b_new->find(mockup_data[0].first) == mockup_data[0].second);
}

TEST_F(AuxParallelGrowingHTGTest, testHTCustodianCopyConstructor) {
    Aux::ParallelHashMap phm_a;
    auto mockup_data = generateMockupData();

    auto handle_a = phm_a.make_handle();

    bool all_inserted_successful = insert_data(*handle_a.get(), mockup_data);

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, handle_a->hashtable()));

    Aux::ParallelHashMap phm_b(phm_a);

    auto handle_b = phm_b.make_handle();

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, handle_b->hashtable()));
}

TEST_F(AuxParallelGrowingHTGTest, testHTCustodianMoveConstructor) {

    Aux::ParallelHashMap phm_a;
    auto mockup_data = generateMockupData();

    auto handle_a = phm_a.make_handle();

    bool all_inserted_successful = insert_data(*handle_a.get(), mockup_data);

    Aux::ParallelHashMap phm_b(std::move(phm_a));

    auto handle_b = phm_b.make_handle();

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, handle_b->hashtable()));
}

TEST_F(AuxParallelGrowingHTGTest, testHTCustodianCopyAssignment) {
    Aux::ParallelHashMap phm_a;
    auto mockup_data = generateMockupData();

    auto handle_a = phm_a.make_handle();

    bool all_inserted_successful = insert_data(*handle_a.get(), mockup_data);

    Aux::ParallelHashMap phm_b;
    phm_b = phm_a;

    auto handle_b = phm_b.make_handle();

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, handle_b->hashtable()));
}

TEST_F(AuxParallelGrowingHTGTest, testHTCustodianSingleThread) {
    Aux::setNumberOfThreads(1);

    Aux::ParallelHashMap phm{32};

    constexpr size_t fill_64bit_values = 10254;
    MockupData mockup_data = generateMockupData(fill_64bit_values);

#pragma omp parallel num_threads(1) shared(phm)
    {
        auto handle = phm.make_handle();
        for (auto it = std::cbegin(mockup_data); it != std::cend(mockup_data); ++it) {
            handle->insert(it->first, it->second);
        }
    }

    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, *phm.current_table()));
}

TEST_F(AuxParallelGrowingHTGTest, testHTCustodianDualThread) {
    Aux::setNumberOfThreads(2);

    Aux::ParallelHashMap phm{32};

    constexpr size_t fill_64bit_values = 10254;
    MockupData mockup_data = generateMockupData(fill_64bit_values);
    std::vector<bool> successful_inserts({true, true});

#pragma omp parallel num_threads(2) shared(phm, mockup_data)
    {
        uint32_t const thread_id = omp_get_thread_num();
        auto handle = phm.make_handle();

        auto local_begin = std::cbegin(mockup_data);
        size_t const elements_per_thread = mockup_data.size() / 2;

        if (thread_id != 0) {
            std::advance(local_begin, thread_id * elements_per_thread);
        }

        auto local_end = local_begin;

        if (thread_id == 1) {
            local_end = std::cend(mockup_data);
        } else {
            std::advance(local_end, (thread_id + 1) * elements_per_thread);
        }

        for (auto it = local_begin; it != local_end && successful_inserts[thread_id]; ++it) {
            successful_inserts[thread_id] = handle->insert(it->first, it->second);
        }
    }

    ASSERT_TRUE(successful_inserts[0]);
    ASSERT_TRUE(successful_inserts[1]);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, *phm.current_table()));
}

TEST_F(AuxParallelGrowingHTGTest, testHTCustodianQuadrupleThread) {
    Aux::setNumberOfThreads(4);

    Aux::ParallelHashMap phm{32};

    constexpr size_t fill_64bit_values = 10254;
    MockupData mockup_data = generateMockupData(fill_64bit_values);
    std::vector<bool> successful_inserts({true, true, true, true});

#pragma omp parallel num_threads(4) shared(phm, mockup_data)
    {
        uint32_t const thread_id = omp_get_thread_num();
        auto handle = phm.make_handle();

        auto local_begin = std::cbegin(mockup_data);
        size_t const elements_per_thread = mockup_data.size() / 4;

        if (thread_id != 0) {
            std::advance(local_begin, thread_id * elements_per_thread);
        }

        auto local_end = local_begin;

        if (thread_id == 3) {
            local_end = std::cend(mockup_data);
        } else {
            std::advance(local_end, elements_per_thread);
        }

        for (auto it = local_begin; it != local_end; ++it) {
            successful_inserts[thread_id] = handle->insert(it->first, it->second);
        }
    }

    ASSERT_TRUE(successful_inserts[0]);
    ASSERT_TRUE(successful_inserts[1]);
    ASSERT_TRUE(successful_inserts[2]);
    ASSERT_TRUE(successful_inserts[3]);
    ASSERT_TRUE(checkEntriesForMockupData(mockup_data, *phm.current_table()));
}

} // namespace NetworKit
