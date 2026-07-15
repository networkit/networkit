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
#include <vector>

#include <networkit/structures/GenericPartition.hpp>

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::SizeIs;
using ::testing::UnorderedElementsAre;

template <typename IndexType>
class GenericPartitionGTest : public ::testing::Test {};

TYPED_TEST_SUITE_P(GenericPartitionGTest);

template <typename IndexType>
GenericPartition<IndexType> makeMergedPartition() {
    GenericPartition<IndexType> partition(IndexType{10});
    partition.allToSingletons();
    partition.mergeSubsets(partition[IndexType{0}], partition[IndexType{9}]);
    partition.mergeSubsets(partition[IndexType{1}], partition[IndexType{8}]);
    partition.mergeSubsets(partition[IndexType{2}], partition[IndexType{7}]);
    partition.mergeSubsets(partition[IndexType{0}], partition[IndexType{1}]);
    partition.mergeSubsets(partition[IndexType{1}], partition[IndexType{2}]);
    return partition;
}

TYPED_TEST_P(GenericPartitionGTest, testSizeConstructorCreatesUnassignedElements) {
    constexpr auto none = GenericPartition<TypeParam>::noneIndex;

    GenericPartition<TypeParam> partition(TypeParam{4});

    EXPECT_EQ(partition.numberOfElements(), TypeParam{4});
    EXPECT_EQ(partition.lowerBound(), TypeParam{0});
    EXPECT_EQ(partition.upperBound(), TypeParam{1});
    EXPECT_THAT(partition.getSubsetIds(), IsEmpty());
    EXPECT_THAT(partition.getVector(), ElementsAre(none, none, none, none));
}

TYPED_TEST_P(GenericPartitionGTest, testSizeAndDefaultValueConstructorInitializesElements) {
    GenericPartition<TypeParam> partitionWithDefaultValue(TypeParam{10}, TypeParam{1});

    EXPECT_EQ(partitionWithDefaultValue.numberOfElements(), TypeParam{10});
    EXPECT_EQ(partitionWithDefaultValue.lowerBound(), TypeParam{0});
    EXPECT_EQ(partitionWithDefaultValue.upperBound(), TypeParam{1});
    EXPECT_THAT(partitionWithDefaultValue.getVector(),
                ElementsAre(TypeParam{1}, TypeParam{1}, TypeParam{1}, TypeParam{1}, TypeParam{1},
                            TypeParam{1}, TypeParam{1}, TypeParam{1}, TypeParam{1}, TypeParam{1}));
}

TYPED_TEST_P(GenericPartitionGTest, testVectorConstructorInitializesFromVector) {
    constexpr auto none = GenericPartition<TypeParam>::noneIndex;

    const std::vector<TypeParam> unassignedElements(3, none);
    GenericPartition<TypeParam> partitionFromVector(unassignedElements);

    EXPECT_EQ(partitionFromVector.numberOfElements(), TypeParam{3});
    EXPECT_EQ(partitionFromVector.lowerBound(), TypeParam{0});
    EXPECT_EQ(partitionFromVector.upperBound(), TypeParam{1});
    EXPECT_THAT(partitionFromVector.getSubsetIds(), IsEmpty());
    EXPECT_THAT(partitionFromVector.getVector(), ElementsAre(none, none, none));
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

TYPED_TEST_P(GenericPartitionGTest, testAddToSubset) {
    constexpr auto none = GenericPartition<TypeParam>::noneIndex;

    GenericPartition<TypeParam> partition(TypeParam{10});

    partition.addToSubset(TypeParam{0}, TypeParam{4});
    partition.addToSubset(TypeParam{0}, TypeParam{7});
    EXPECT_THAT(partition[TypeParam{4}], Eq(partition[TypeParam{7}]));
    EXPECT_THAT(partition[TypeParam{1}], Eq(none));
    EXPECT_THAT(partition[TypeParam{2}], Eq(none));

    partition.setUpperBound(TypeParam{2});
    partition.addToSubset(TypeParam{1}, TypeParam{1});
    partition.addToSubset(TypeParam{1}, TypeParam{2});
    EXPECT_THAT(partition[TypeParam{4}], Eq(partition[TypeParam{7}]));
    EXPECT_THAT(partition[TypeParam{1}], Eq(partition[TypeParam{2}]));
    EXPECT_THAT(partition[TypeParam{1}], Eq(TypeParam{1}));
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

TYPED_TEST_P(GenericPartitionGTest, testNumberOfSubsetsAfterSingletons) {
    constexpr auto elements = TypeParam{517};

    GenericPartition<TypeParam> partition(elements);
    partition.allToSingletons();
    EXPECT_THAT(partition.numberOfSubsets(), Eq(count{elements}));
}

TYPED_TEST_P(GenericPartitionGTest, testNumberOfSubsetsAfterMovingAndMerging) {
    GenericPartition<TypeParam> oneSubsetPartition(TypeParam{10});
    oneSubsetPartition.allToSingletons();

    for (TypeParam element = 0; element < oneSubsetPartition.numberOfElements(); ++element) {
        oneSubsetPartition.moveToSubset(oneSubsetPartition[TypeParam{0}], element);
    }

    EXPECT_THAT(oneSubsetPartition.numberOfSubsets(), Eq(count{1}));

    GenericPartition<TypeParam> pairwiseMergedPartition(TypeParam{10});
    pairwiseMergedPartition.allToSingletons();
    for (TypeParam element = 0; element < pairwiseMergedPartition.numberOfElements();
         element += TypeParam{2}) {
        pairwiseMergedPartition.mergeSubsets(pairwiseMergedPartition[element],
                                             pairwiseMergedPartition[element + TypeParam{1}]);
    }

    EXPECT_THAT(pairwiseMergedPartition.numberOfSubsets(), Eq(count{5}));
    EXPECT_THAT(pairwiseMergedPartition.upperBound(), Eq(TypeParam{15}));
}

TYPED_TEST_P(GenericPartitionGTest, testBoundsAfterSingletons) {
    GenericPartition<TypeParam> partition(TypeParam{10});
    partition.allToSingletons();

    EXPECT_THAT(partition.lowerBound(), Eq(TypeParam{0}));
    EXPECT_THAT(partition.upperBound(), Eq(TypeParam{10}));
}

TYPED_TEST_P(GenericPartitionGTest, testContainsAssignedElementsOnly) {
    GenericPartition<TypeParam> unassignedPartition(TypeParam{10});
    EXPECT_THAT(unassignedPartition.contains(TypeParam{0}), Eq(false));
    EXPECT_THAT(unassignedPartition.contains(TypeParam{9}), Eq(false));

    unassignedPartition.allToSingletons();
    EXPECT_THAT(unassignedPartition.contains(TypeParam{8}), Eq(true));
    EXPECT_THAT(unassignedPartition.contains(TypeParam{15}), Eq(false));

    auto mergedPartition = makeMergedPartition<TypeParam>();
    mergedPartition.toSingleton(TypeParam{9});
    EXPECT_THAT(mergedPartition.contains(TypeParam{8}), Eq(true));
    EXPECT_THAT(mergedPartition.contains(TypeParam{9}), Eq(true));
}

TYPED_TEST_P(GenericPartitionGTest, testMergedPartitionStructure) {
    auto partition = makeMergedPartition<TypeParam>();

    EXPECT_THAT(partition.numberOfElements(), Eq(TypeParam{10}));
    EXPECT_THAT(partition.inSameSubset(TypeParam{9}, TypeParam{7}), Eq(true));
    EXPECT_THAT(partition.inSameSubset(TypeParam{9}, TypeParam{4}), Eq(false));
    EXPECT_THAT(partition.getMembers(partition[TypeParam{0}]),
                ElementsAre(TypeParam{0}, TypeParam{1}, TypeParam{2}, TypeParam{7}, TypeParam{8},
                            TypeParam{9}));
    EXPECT_THAT(partition.subsetSizes(),
                ElementsAre(count{1}, count{1}, count{1}, count{1}, count{6}));
    EXPECT_THAT(partition.subsetSizeMap(),
                ElementsAre(Pair(TypeParam{3}, count{1}), Pair(TypeParam{4}, count{1}),
                            Pair(TypeParam{5}, count{1}), Pair(TypeParam{6}, count{1}),
                            Pair(partition[TypeParam{0}], count{6})));
}

TYPED_TEST_P(GenericPartitionGTest, testCompactPreservesMergedPartitionStructure) {
    auto partition = makeMergedPartition<TypeParam>();

    partition.compact();

    EXPECT_THAT(partition.upperBound(), Eq(TypeParam{5}));
    EXPECT_THAT(partition.getSubsets(),
                UnorderedElementsAre(ElementsAre(TypeParam{0}, TypeParam{1}, TypeParam{2},
                                                 TypeParam{7}, TypeParam{8}, TypeParam{9}),
                                     ElementsAre(TypeParam{3}), ElementsAre(TypeParam{4}),
                                     ElementsAre(TypeParam{5}), ElementsAre(TypeParam{6})));
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

REGISTER_TYPED_TEST_SUITE_P(GenericPartitionGTest, testSizeConstructorCreatesUnassignedElements,
                            testSizeAndDefaultValueConstructorInitializesElements,
                            testVectorConstructorInitializesFromVector,
                            testSingletonsUseTemplateIndexType, testMoveMergeAndMembers,
                            testCompactRenumbersSparseSubsetIds,
                            testRemoveAndResetUseTypedNoneSentinel, testAddToSubset,
                            testNumberOfSubsetsAfterSingletons,
                            testNumberOfSubsetsAfterMovingAndMerging, testBoundsAfterSingletons,
                            testContainsAssignedElementsOnly, testMergedPartitionStructure,
                            testCompactPreservesMergedPartitionStructure);

using GenericPartitionTestTypes =
    ::testing::Types<index, std::uint32_t, std::uint16_t, std::int64_t, std::int32_t, std::int16_t>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestGenericPartition, GenericPartitionGTest,
                               GenericPartitionTestTypes, );

static_assert(std::same_as<Partition, GenericPartition<index>>);

} // namespace
} // namespace NetworKit
