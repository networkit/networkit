/*
 * GenericPartitionGTest.cpp
 *
 *  Created on: 10.07.2026
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <set>
#include <type_traits>
#include <vector>

#include <networkit/structures/GenericPartition.hpp>

namespace NetworKit {
namespace {

template <typename IndexT>
struct GenericPartitionConfig {
    using IndexType = IndexT;
};

template <typename TestT>
class GenericPartitionGTest : public ::testing::Test {
public:
    using IndexType = typename TestT::IndexType;
};

TYPED_TEST_SUITE_P(GenericPartitionGTest);

TYPED_TEST_P(GenericPartitionGTest, testConstructorCreatesUnassignedElements) {
    using ::testing::ElementsAre;
    using ::testing::IsEmpty;

    using IndexType = typename TestFixture::IndexType;
    constexpr auto none = std::numeric_limits<IndexType>::max();

    GenericPartition<IndexType> partition(IndexType{4});

    EXPECT_EQ(partition.numberOfElements(), IndexType{4});
    EXPECT_EQ(partition.lowerBound(), IndexType{0});
    EXPECT_EQ(partition.upperBound(), IndexType{1});
    EXPECT_THAT(partition.getSubsetIds(), IsEmpty());
    EXPECT_THAT(partition.getVector(), ElementsAre(none, none, none, none));
}

TYPED_TEST_P(GenericPartitionGTest, testSingletonsUseTemplateIndexType) {
    using ::testing::ElementsAre;
    using ::testing::SizeIs;

    using IndexType = typename TestFixture::IndexType;

    GenericPartition<IndexType> partition(IndexType{5});
    partition.allToSingletons();

    EXPECT_EQ(partition.numberOfElements(), IndexType{5});
    EXPECT_THAT(partition.getSubsetIds(),
                ElementsAre(IndexType{0}, IndexType{1}, IndexType{2}, IndexType{3}, IndexType{4}));
    EXPECT_THAT(partition.getSubsets(), SizeIs(5));

    for (IndexType i = 0; i < partition.numberOfElements(); ++i) {
        EXPECT_EQ(partition.subsetOf(i), i);
        EXPECT_TRUE(partition.contains(i));
    }
}

TYPED_TEST_P(GenericPartitionGTest, testMoveMergeAndMembers) {
    using ::testing::ElementsAre;
    using ::testing::UnorderedElementsAre;

    using IndexType = typename TestFixture::IndexType;

    GenericPartition<IndexType> partition(IndexType{6});
    partition.allToSingletons();

    partition.moveToSubset(IndexType{1}, IndexType{4});
    EXPECT_TRUE(partition.inSameSubset(IndexType{1}, IndexType{4}));
    EXPECT_THAT(partition.getMembers(IndexType{1}), ElementsAre(IndexType{1}, IndexType{4}));

    const auto merged = partition.mergeSubsets(IndexType{1}, IndexType{3});
    EXPECT_EQ(partition.subsetOf(IndexType{1}), merged);
    EXPECT_EQ(partition.subsetOf(IndexType{3}), merged);
    EXPECT_EQ(partition.subsetOf(IndexType{4}), merged);
    EXPECT_THAT(partition.getMembers(merged),
                UnorderedElementsAre(IndexType{1}, IndexType{3}, IndexType{4}));
}

TYPED_TEST_P(GenericPartitionGTest, testCompactRenumbersSparseSubsetIds) {
    using ::testing::ElementsAre;
    using ::testing::Pair;

    using IndexType = typename TestFixture::IndexType;
    constexpr auto none = std::numeric_limits<IndexType>::max();

    GenericPartition<IndexType> partition(IndexType{5}, none);
    partition.setUpperBound(IndexType{21});
    partition.addToSubset(IndexType{20}, IndexType{0});
    partition.addToSubset(IndexType{7}, IndexType{1});
    partition.addToSubset(IndexType{20}, IndexType{2});
    partition.addToSubset(IndexType{7}, IndexType{3});

    partition.compact();

    EXPECT_EQ(partition.numberOfSubsets(), count{2});
    EXPECT_EQ(partition.upperBound(), IndexType{2});
    EXPECT_THAT(partition.subsetSizeMap(),
                ElementsAre(Pair(IndexType{0}, count{2}), Pair(IndexType{1}, count{2})));
    EXPECT_EQ(partition.subsetOf(IndexType{4}), none);
}

TYPED_TEST_P(GenericPartitionGTest, testRemoveAndResetUseTypedNoneSentinel) {
    using ::testing::ElementsAre;
    using ::testing::IsEmpty;

    using IndexType = typename TestFixture::IndexType;
    constexpr auto none = std::numeric_limits<IndexType>::max();

    GenericPartition<IndexType> partition(IndexType{3});
    partition.allToOnePartition();
    partition.remove(IndexType{1});

    EXPECT_FALSE(partition.contains(IndexType{1}));
    EXPECT_EQ(partition.subsetOf(IndexType{1}), none);

    partition.reset(IndexType{2}, none);
    EXPECT_THAT(partition.getSubsetIds(), IsEmpty());
    EXPECT_THAT(partition.getVector(), ElementsAre(none, none));
}

REGISTER_TYPED_TEST_SUITE_P(GenericPartitionGTest, testConstructorCreatesUnassignedElements,
                            testSingletonsUseTemplateIndexType, testMoveMergeAndMembers,
                            testCompactRenumbersSparseSubsetIds,
                            testRemoveAndResetUseTypedNoneSentinel);

using GenericPartitionTestTypes =
    ::testing::Types<GenericPartitionConfig<index>,
                     GenericPartitionConfig<std::uint32_t>,
                     GenericPartitionConfig<std::uint16_t>,
                     GenericPartitionConfig<std::int64_t>,
                     GenericPartitionConfig<std::int32_t>,
                     GenericPartitionConfig<std::int16_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestGenericPartition, GenericPartitionGTest,
                               GenericPartitionTestTypes, );

static_assert(std::same_as<Partition, GenericPartition<index>>);

} // namespace
} // namespace NetworKit
