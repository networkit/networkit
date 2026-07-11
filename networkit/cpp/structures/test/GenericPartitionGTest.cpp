/*
 * GenericPartitionGTest.cpp
 *
 *  Created on: 10.07.2026
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <concepts>
#include <cstdint>

#include <networkit/structures/GenericPartition.hpp>

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::SizeIs;
using ::testing::UnorderedElementsAre;

template <typename IndexType>
class GenericPartitionGTest : public ::testing::Test {};

TYPED_TEST_SUITE_P(GenericPartitionGTest);

TYPED_TEST_P(GenericPartitionGTest, testConstructorCreatesUnassignedElements) {
    constexpr auto none = GenericPartition<TypeParam>::noneIndex;

    GenericPartition<TypeParam> partition(TypeParam{4});

    EXPECT_EQ(partition.numberOfElements(), TypeParam{4});
    EXPECT_EQ(partition.lowerBound(), TypeParam{0});
    EXPECT_EQ(partition.upperBound(), TypeParam{1});
    EXPECT_THAT(partition.getSubsetIds(), IsEmpty());
    EXPECT_THAT(partition.getVector(), ElementsAre(none, none, none, none));
}

TYPED_TEST_P(GenericPartitionGTest, testSingletonsUseTemplateIndexType) {
    GenericPartition<TypeParam> partition(TypeParam{5});
    partition.allToSingletons();

    EXPECT_EQ(partition.numberOfElements(), TypeParam{5});
    EXPECT_THAT(partition.getSubsetIds(),
                ElementsAre(TypeParam{0}, TypeParam{1}, TypeParam{2}, TypeParam{3}, TypeParam{4}));
    EXPECT_THAT(partition.getSubsets(), SizeIs(5));

    for (TypeParam i = 0; i < partition.numberOfElements(); ++i) {
        EXPECT_EQ(partition.subsetOf(i), i);
        EXPECT_TRUE(partition.contains(i));
    }
}

TYPED_TEST_P(GenericPartitionGTest, testMoveMergeAndMembers) {
    GenericPartition<TypeParam> partition(TypeParam{6});
    partition.allToSingletons();

    partition.moveToSubset(TypeParam{1}, TypeParam{4});
    EXPECT_TRUE(partition.inSameSubset(TypeParam{1}, TypeParam{4}));
    EXPECT_THAT(partition.getMembers(TypeParam{1}), ElementsAre(TypeParam{1}, TypeParam{4}));

    const auto merged = partition.mergeSubsets(TypeParam{1}, TypeParam{3});
    EXPECT_EQ(partition.subsetOf(TypeParam{1}), merged);
    EXPECT_EQ(partition.subsetOf(TypeParam{3}), merged);
    EXPECT_EQ(partition.subsetOf(TypeParam{4}), merged);
    EXPECT_THAT(partition.getMembers(merged),
                UnorderedElementsAre(TypeParam{1}, TypeParam{3}, TypeParam{4}));
}

TYPED_TEST_P(GenericPartitionGTest, testCompactRenumbersSparseSubsetIds) {
    constexpr auto none = GenericPartition<TypeParam>::noneIndex;

    GenericPartition<TypeParam> partition(TypeParam{5}, none);
    partition.setUpperBound(TypeParam{21});
    partition.addToSubset(TypeParam{20}, TypeParam{0});
    partition.addToSubset(TypeParam{7}, TypeParam{1});
    partition.addToSubset(TypeParam{20}, TypeParam{2});
    partition.addToSubset(TypeParam{7}, TypeParam{3});

    partition.compact();

    EXPECT_EQ(partition.numberOfSubsets(), count{2});
    EXPECT_EQ(partition.upperBound(), TypeParam{2});
    EXPECT_THAT(partition.subsetSizeMap(),
                ElementsAre(Pair(TypeParam{0}, count{2}), Pair(TypeParam{1}, count{2})));
    EXPECT_EQ(partition.subsetOf(TypeParam{4}), none);
}

TYPED_TEST_P(GenericPartitionGTest, testRemoveAndResetUseTypedNoneSentinel) {
    constexpr auto none = GenericPartition<TypeParam>::noneIndex;

    GenericPartition<TypeParam> partition(TypeParam{3});
    partition.allToOnePartition();
    partition.remove(TypeParam{1});

    EXPECT_FALSE(partition.contains(TypeParam{1}));
    EXPECT_EQ(partition.subsetOf(TypeParam{1}), none);

    partition.reset(TypeParam{2}, none);
    EXPECT_THAT(partition.getSubsetIds(), IsEmpty());
    EXPECT_THAT(partition.getVector(), ElementsAre(none, none));
}

REGISTER_TYPED_TEST_SUITE_P(GenericPartitionGTest, testConstructorCreatesUnassignedElements,
                            testSingletonsUseTemplateIndexType, testMoveMergeAndMembers,
                            testCompactRenumbersSparseSubsetIds,
                            testRemoveAndResetUseTypedNoneSentinel);

using GenericPartitionTestTypes =
    ::testing::Types<index, std::uint32_t, std::uint16_t, std::int64_t, std::int32_t, std::int16_t>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestGenericPartition, GenericPartitionGTest,
                               GenericPartitionTestTypes, );

static_assert(std::same_as<Partition, GenericPartition<index>>);

} // namespace
} // namespace NetworKit
