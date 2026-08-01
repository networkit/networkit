/*
 * GenericCoverGTest.cpp
 *
 *  Created on: 12.12.2013
 *      Author: Maximilian Vogel (uocvf@student.kit.edu)
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <concepts>
#include <cstdint>
#include <set>

#include <networkit/structures/GenericCover.hpp>

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::Pair;

template <typename IndexType>
class GenericCoverGTest : public ::testing::Test {};

TYPED_TEST_SUITE_P(GenericCoverGTest);

template <typename IndexType>
GenericCover<IndexType> makeOverlappingCover() {
    constexpr IndexType n{10};
    GenericCover<IndexType> cover(n);
    for (IndexType i = 0; i < n; i += IndexType{2}) {
        IndexType sid = cover.toSingleton(i);
        cover.toSingleton(i + IndexType{1});
        cover.addToSubset(sid, i + IndexType{1});
    }
    for (IndexType i = 0; i < n; i++) {
        cover.addToSubset(i + IndexType{1}, IndexType{0});
    }
    return cover;
}

TYPED_TEST_P(GenericCoverGTest, testConstructor) {
    GenericCover<TypeParam> cover(TypeParam{10});
    EXPECT_EQ(cover.lowerBound(), TypeParam{0});
    EXPECT_EQ(cover.upperBound(), TypeParam{1});
}

TYPED_TEST_P(GenericCoverGTest, testPartitionConstructor) {
    GenericPartition<TypeParam> partition(TypeParam{4});
    partition.allToSingletons();

    GenericCover<TypeParam> cover(partition);

    EXPECT_EQ(cover.numberOfElements(), count{4});
    EXPECT_EQ(cover.upperBound(), TypeParam{4});
    EXPECT_THAT(cover.subsetsOf(TypeParam{0}), ElementsAre(TypeParam{0}));
    EXPECT_THAT(cover.subsetsOf(TypeParam{3}), ElementsAre(TypeParam{3}));
}

TYPED_TEST_P(GenericCoverGTest, testAllToSingletonsAndUpperBound) {
    GenericCover<TypeParam> cover(TypeParam{10});
    EXPECT_EQ(cover.upperBound(), TypeParam{1});
    cover.allToSingletons();
    EXPECT_EQ(cover.lowerBound(), TypeParam{0});
    EXPECT_EQ(cover.upperBound(), TypeParam{11});
}

TYPED_TEST_P(GenericCoverGTest, testContains) {
    GenericCover<TypeParam> cover(TypeParam{10});
    cover.toSingleton(TypeParam{0});
    EXPECT_TRUE(cover.contains(TypeParam{0}));
    EXPECT_FALSE(cover.contains(TypeParam{1}));
}

TYPED_TEST_P(GenericCoverGTest, testUpperBoundAfterMerges) {
    GenericCover<TypeParam> cover = makeOverlappingCover<TypeParam>();
    cover.mergeSubsets(TypeParam{1}, TypeParam{3});
    cover.mergeSubsets(TypeParam{5}, TypeParam{11});
    EXPECT_EQ(cover.upperBound(), TypeParam{13});
}

TYPED_TEST_P(GenericCoverGTest, testToSingleton) {
    GenericCover<TypeParam> cover(TypeParam{10});
    cover.allToSingletons();

    EXPECT_THAT(cover[TypeParam{0}], ElementsAre(TypeParam{1}));
    cover.addToSubset(TypeParam{5}, TypeParam{0});
    cover.addToSubset(TypeParam{2}, TypeParam{0});
    cover.addToSubset(TypeParam{3}, TypeParam{0});
    cover.addToSubset(TypeParam{4}, TypeParam{0});
    cover.addToSubset(TypeParam{0}, TypeParam{1});
    cover.toSingleton(TypeParam{0});

    EXPECT_THAT(cover[TypeParam{0}], ElementsAre(TypeParam{11}));
}

TYPED_TEST_P(GenericCoverGTest, testAddToSubset) {
    GenericCover<TypeParam> cover(TypeParam{10});
    cover.addToSubset(TypeParam{0}, TypeParam{0});
    cover.addToSubset(TypeParam{0}, TypeParam{1});

    EXPECT_TRUE(cover.inSameSubset(TypeParam{0}, TypeParam{1}));
    EXPECT_THAT(cover[TypeParam{0}], ElementsAre(TypeParam{0}));
    EXPECT_THAT(cover[TypeParam{1}], ElementsAre(TypeParam{0}));
}

TYPED_TEST_P(GenericCoverGTest, testAddToSubset2) {
    GenericCover<TypeParam> cover(TypeParam{10});
    TypeParam sid = cover.toSingleton(TypeParam{0});
    cover.addToSubset(sid, TypeParam{5});
    EXPECT_TRUE(cover.inSameSubset(TypeParam{0}, TypeParam{5}));
}

TYPED_TEST_P(GenericCoverGTest, testMoveToSubset) {
    GenericCover<TypeParam> cover(TypeParam{10});
    cover.allToSingletons();
    cover.addToSubset(TypeParam{5}, TypeParam{0});
    cover.addToSubset(TypeParam{2}, TypeParam{0});
    cover.addToSubset(TypeParam{3}, TypeParam{0});
    cover.addToSubset(TypeParam{4}, TypeParam{0});
    cover.addToSubset(TypeParam{0}, TypeParam{1});
    cover.moveToSubset(TypeParam{8}, TypeParam{0});
    EXPECT_THAT(cover[TypeParam{0}], ElementsAre(TypeParam{8}));
}

TYPED_TEST_P(GenericCoverGTest, testSubsetSizesWithUnassignedElements) {
    constexpr TypeParam n{10};
    GenericCover<TypeParam> cover(n);
    for (TypeParam i = 0; i < n; i += TypeParam{2}) {
        cover.toSingleton(i);
    }
    EXPECT_THAT(cover.subsetSizes(), ElementsAre(count{1}, count{1}, count{1}, count{1}, count{1}));
}

TYPED_TEST_P(GenericCoverGTest, testSubsetSizesTrivial) {
    constexpr TypeParam n{10};
    GenericCover<TypeParam> cover(n);
    for (TypeParam i = 0; i < n; i++) {
        cover.toSingleton(i);
    }
    EXPECT_THAT(cover.subsetSizes(), ElementsAre(count{1}, count{1}, count{1}, count{1}, count{1},
                                                 count{1}, count{1}, count{1}, count{1}, count{1}));
}

TYPED_TEST_P(GenericCoverGTest, testSubsetSizesTrivial2) {
    constexpr TypeParam n{10};
    GenericCover<TypeParam> cover(n);
    for (TypeParam i = 0; i < n; i += TypeParam{2}) {
        cover.toSingleton(i);
    }
    for (TypeParam i = 1; i < n; i += TypeParam{2}) {
        cover.addToSubset(i / TypeParam{2} + TypeParam{1}, i);
    }
    EXPECT_THAT(cover.subsetSizes(), ElementsAre(count{2}, count{2}, count{2}, count{2}, count{2}));
}

TYPED_TEST_P(GenericCoverGTest, testSubsetSizesAssignedToMultipleSubsets) {
    constexpr TypeParam n{10};
    GenericCover<TypeParam> cover(n);
    for (TypeParam i = 0; i < n; i++) {
        cover.toSingleton(i);
    }
    for (TypeParam i = 1; i < n; i += TypeParam{2}) {
        cover.addToSubset(i, i);
    }
    EXPECT_THAT(cover.subsetSizes(), ElementsAre(count{2}, count{1}, count{2}, count{1}, count{2},
                                                 count{1}, count{2}, count{1}, count{2}, count{1}));
}

TYPED_TEST_P(GenericCoverGTest, testSubsetSizesAssignedToMultipleSubsets2) {
    GenericCover<TypeParam> cover = makeOverlappingCover<TypeParam>();

    EXPECT_THAT(cover.subsetSizes(), ElementsAre(count{2}, count{2}, count{3}, count{2}, count{3},
                                                 count{2}, count{3}, count{2}, count{3}, count{2}));
}

TYPED_TEST_P(GenericCoverGTest, testSubsetSizeMapMultipleSets) {
    constexpr TypeParam n{10};
    GenericCover<TypeParam> cover(n);
    for (TypeParam i = 0; i < n; i++) {
        cover.toSingleton(i);
    }
    for (TypeParam i = 1; i < n; i += TypeParam{2}) {
        cover.addToSubset(i, i);
    }

    EXPECT_THAT(cover.subsetSizeMap(),
                ElementsAre(Pair(TypeParam{1}, count{2}), Pair(TypeParam{2}, count{1}),
                            Pair(TypeParam{3}, count{2}), Pair(TypeParam{4}, count{1}),
                            Pair(TypeParam{5}, count{2}), Pair(TypeParam{6}, count{1}),
                            Pair(TypeParam{7}, count{2}), Pair(TypeParam{8}, count{1}),
                            Pair(TypeParam{9}, count{2}), Pair(TypeParam{10}, count{1})));
}

TYPED_TEST_P(GenericCoverGTest, testMergeSubsetsAndGetMembers) {
    GenericCover<TypeParam> cover = makeOverlappingCover<TypeParam>();

    cover.mergeSubsets(TypeParam{1}, TypeParam{3});
    cover.mergeSubsets(TypeParam{5}, TypeParam{11});

    // remaining subset IDs 2,4,6,7,8,9,10,12
    // their sizes          2,2,2,3,2,3,2,6
    EXPECT_THAT(cover.subsetSizes(), ElementsAre(count{2}, count{2}, count{2}, count{3}, count{2},
                                                 count{3}, count{2}, count{6}));
    EXPECT_THAT(cover.getMembers(TypeParam{12}),
                ElementsAre(TypeParam{0}, TypeParam{1}, TypeParam{2}, TypeParam{3}, TypeParam{4},
                            TypeParam{5}));
}

TYPED_TEST_P(GenericCoverGTest, testNumberOfSubsets) {
    GenericCover<TypeParam> cover = makeOverlappingCover<TypeParam>();

    EXPECT_EQ(cover.numberOfSubsets(), count{10});
    cover.mergeSubsets(TypeParam{1}, TypeParam{2});
    cover.mergeSubsets(TypeParam{3}, TypeParam{11});
    EXPECT_EQ(cover.numberOfSubsets(), count{8});
}

TYPED_TEST_P(GenericCoverGTest, testSubsetsOf) {
    GenericCover<TypeParam> cover = makeOverlappingCover<TypeParam>();

    EXPECT_THAT(cover.subsetsOf(TypeParam{0}),
                ElementsAre(TypeParam{1}, TypeParam{2}, TypeParam{3}, TypeParam{4}, TypeParam{5},
                            TypeParam{6}, TypeParam{7}, TypeParam{8}, TypeParam{9}, TypeParam{10}));
    EXPECT_THAT(cover.subsetsOf(TypeParam{3}), ElementsAre(TypeParam{3}, TypeParam{4}));
    cover.mergeSubsets(TypeParam{1}, TypeParam{3});
    cover.mergeSubsets(TypeParam{5}, TypeParam{11});
    EXPECT_THAT(cover.subsetsOf(TypeParam{0}),
                ElementsAre(TypeParam{2}, TypeParam{4}, TypeParam{6}, TypeParam{7}, TypeParam{8},
                            TypeParam{9}, TypeParam{10}, TypeParam{12}));
}

TYPED_TEST_P(GenericCoverGTest, testInSameSubset) {
    constexpr TypeParam n{10};
    GenericCover<TypeParam> cover(n);
    for (TypeParam i = 0; i < n; i += TypeParam{2}) {
        TypeParam sid = cover.toSingleton(i);
        cover.toSingleton(i + TypeParam{1});
        cover.addToSubset(sid, i + TypeParam{1});
    }
    EXPECT_TRUE(cover.inSameSubset(TypeParam{0}, TypeParam{1}));
    EXPECT_FALSE(cover.inSameSubset(TypeParam{0}, TypeParam{2}));
    EXPECT_FALSE(cover.inSameSubset(TypeParam{1}, TypeParam{5}));
    cover.mergeSubsets(TypeParam{1}, TypeParam{3});
    EXPECT_TRUE(cover.inSameSubset(TypeParam{0}, TypeParam{1}));
    EXPECT_TRUE(cover.inSameSubset(TypeParam{0}, TypeParam{2}));
    EXPECT_FALSE(cover.inSameSubset(TypeParam{1}, TypeParam{5}));
    cover.mergeSubsets(TypeParam{5}, TypeParam{11});
    EXPECT_TRUE(cover.inSameSubset(TypeParam{0}, TypeParam{1}));
    EXPECT_TRUE(cover.inSameSubset(TypeParam{0}, TypeParam{2}));
    EXPECT_TRUE(cover.inSameSubset(TypeParam{1}, TypeParam{5}));
}

REGISTER_TYPED_TEST_SUITE_P(GenericCoverGTest, testConstructor, testPartitionConstructor,
                            testAllToSingletonsAndUpperBound, testContains,
                            testUpperBoundAfterMerges, testToSingleton, testAddToSubset,
                            testAddToSubset2, testMoveToSubset,
                            testSubsetSizesWithUnassignedElements, testSubsetSizesTrivial,
                            testSubsetSizesTrivial2, testSubsetSizesAssignedToMultipleSubsets,
                            testSubsetSizesAssignedToMultipleSubsets2,
                            testSubsetSizeMapMultipleSets, testMergeSubsetsAndGetMembers,
                            testNumberOfSubsets, testSubsetsOf, testInSameSubset);

using GenericCoverTestTypes =
    ::testing::Types<index, std::uint32_t, std::uint16_t, std::int64_t, std::int32_t, std::int16_t>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestGenericCover, GenericCoverGTest, GenericCoverTestTypes, );

static_assert(std::same_as<Cover, GenericCover<index>>);

} // namespace
} // namespace NetworKit
