/*
 * GenericUnionFindGTest.cpp
 *
 *  Created on: 04.12.2013
 *      Author: Maximilian Vogel (uocvf@student.kit.edu)
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <concepts>
#include <cstdint>

#include <networkit/structures/GenericUnionFind.hpp>

namespace NetworKit {
namespace {

using ::testing::Eq;
using ::testing::Ne;

template <typename IndexType>
class GenericUnionFindGTest : public testing::Test {};

TYPED_TEST_SUITE_P(GenericUnionFindGTest);

TYPED_TEST_P(GenericUnionFindGTest, testAllToSingletons) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{10});
    unionFind.allToSingletons();

    for (TypeParam i = 0; i < TypeParam{10}; ++i) {
        EXPECT_THAT(unionFind.find(i),
                    Ne(unionFind.find(static_cast<TypeParam>((i + TypeParam{1}) % TypeParam{10}))));
    }
}

TYPED_TEST_P(GenericUnionFindGTest, testEmptyUnionFindConvertsToEmptyPartition) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{0});

    const auto partition = unionFind.toPartition();

    EXPECT_THAT(partition.numberOfElements(), Eq(TypeParam{0}));
    EXPECT_THAT(partition.numberOfSubsets(), Eq(count{0}));
}

TYPED_TEST_P(GenericUnionFindGTest, testSingleElementAndSelfMerge) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{1});

    unionFind.merge(TypeParam{0}, TypeParam{0});

    EXPECT_THAT(unionFind.find(TypeParam{0}), Eq(TypeParam{0}));

    const auto partition = unionFind.toPartition();
    EXPECT_THAT(partition.numberOfElements(), Eq(TypeParam{1}));
    EXPECT_THAT(partition.numberOfSubsets(), Eq(count{1}));
    EXPECT_THAT(partition.subsetOf(TypeParam{0}), Eq(TypeParam{0}));
}

TYPED_TEST_P(GenericUnionFindGTest, testMergeSimple) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{10});

    unionFind.merge(TypeParam{3}, TypeParam{8});

    EXPECT_THAT(unionFind.find(TypeParam{8}), Eq(unionFind.find(TypeParam{3})));
}

TYPED_TEST_P(GenericUnionFindGTest, testRepeatedMergeOfSameElementsDoesNotChangeSets) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{4});

    unionFind.merge(TypeParam{0}, TypeParam{1});
    const auto representative = unionFind.find(TypeParam{0});

    unionFind.merge(TypeParam{0}, TypeParam{1});
    unionFind.merge(TypeParam{1}, TypeParam{0});

    EXPECT_THAT(unionFind.find(TypeParam{0}), Eq(representative));
    EXPECT_THAT(unionFind.find(TypeParam{1}), Eq(representative));
    EXPECT_THAT(unionFind.find(TypeParam{2}), Ne(representative));
    EXPECT_THAT(unionFind.find(TypeParam{3}), Ne(representative));
}

TYPED_TEST_P(GenericUnionFindGTest, testMergeSubsets) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{10});
    unionFind.allToSingletons();

    unionFind.merge(TypeParam{0}, TypeParam{9});
    unionFind.merge(TypeParam{1}, TypeParam{8});
    unionFind.merge(TypeParam{2}, TypeParam{7});
    unionFind.merge(TypeParam{0}, TypeParam{1});
    unionFind.merge(TypeParam{1}, TypeParam{2});

    EXPECT_THAT(unionFind.find(TypeParam{9}), Eq(unionFind.find(TypeParam{7})));
}

TYPED_TEST_P(GenericUnionFindGTest, testMergeCircular) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{16});

    unionFind.merge(TypeParam{0}, TypeParam{4});
    unionFind.merge(TypeParam{1}, TypeParam{5});
    unionFind.merge(TypeParam{2}, TypeParam{6});
    unionFind.merge(TypeParam{3}, TypeParam{7});
    unionFind.merge(TypeParam{8}, TypeParam{12});
    unionFind.merge(TypeParam{9}, TypeParam{13});
    unionFind.merge(TypeParam{10}, TypeParam{14});
    unionFind.merge(TypeParam{11}, TypeParam{15});

    unionFind.merge(TypeParam{0}, TypeParam{8});
    unionFind.merge(TypeParam{1}, TypeParam{9});
    unionFind.merge(TypeParam{2}, TypeParam{10});
    unionFind.merge(TypeParam{3}, TypeParam{11});
    unionFind.merge(TypeParam{4}, TypeParam{12});
    unionFind.merge(TypeParam{5}, TypeParam{13});
    unionFind.merge(TypeParam{6}, TypeParam{14});
    unionFind.merge(TypeParam{7}, TypeParam{15});

    for (TypeParam i = 0; i < TypeParam{15}; ++i) {
        EXPECT_THAT(unionFind.find(i),
                    Eq(unionFind.find(static_cast<TypeParam>((i + TypeParam{4}) % TypeParam{16}))));
        EXPECT_THAT(unionFind.find(i),
                    Eq(unionFind.find(static_cast<TypeParam>((i + TypeParam{8}) % TypeParam{16}))));
        EXPECT_THAT(unionFind.find(i), Eq(unionFind.find(static_cast<TypeParam>((i + TypeParam{12})
                                                                                % TypeParam{16}))));
        EXPECT_THAT(unionFind.find(i),
                    Ne(unionFind.find(static_cast<TypeParam>((i + TypeParam{1}) % TypeParam{16}))));
        EXPECT_THAT(unionFind.find(i),
                    Ne(unionFind.find(static_cast<TypeParam>((i + TypeParam{2}) % TypeParam{16}))));
        EXPECT_THAT(unionFind.find(i),
                    Ne(unionFind.find(static_cast<TypeParam>((i + TypeParam{3}) % TypeParam{16}))));
    }
}

TYPED_TEST_P(GenericUnionFindGTest, testAllToSingletonsAfterMergesRestoresDistinctSets) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{6});
    unionFind.merge(TypeParam{0}, TypeParam{5});
    unionFind.merge(TypeParam{1}, TypeParam{4});
    unionFind.merge(TypeParam{0}, TypeParam{1});

    unionFind.allToSingletons();

    for (TypeParam i = 0; i < TypeParam{6}; ++i) {
        EXPECT_THAT(unionFind.find(i), Eq(i));
    }
}

TYPED_TEST_P(GenericUnionFindGTest, testToPartitionUsesTemplateIndexType) {
    GenericUnionFind<TypeParam> unionFind(TypeParam{6});
    unionFind.merge(TypeParam{0}, TypeParam{5});
    unionFind.merge(TypeParam{1}, TypeParam{4});
    unionFind.merge(TypeParam{0}, TypeParam{1});

    const auto partition = unionFind.toPartition();

    EXPECT_THAT(partition.numberOfElements(), Eq(TypeParam{6}));
    EXPECT_THAT(partition.numberOfSubsets(), Eq(count{3}));
    EXPECT_THAT(partition.inSameSubset(TypeParam{0}, TypeParam{4}), Eq(true));
    EXPECT_THAT(partition.inSameSubset(TypeParam{5}, TypeParam{1}), Eq(true));
    EXPECT_THAT(partition.inSameSubset(TypeParam{0}, TypeParam{2}), Eq(false));
}

REGISTER_TYPED_TEST_SUITE_P(GenericUnionFindGTest, testAllToSingletons,
                            testEmptyUnionFindConvertsToEmptyPartition,
                            testSingleElementAndSelfMerge, testMergeSimple,
                            testRepeatedMergeOfSameElementsDoesNotChangeSets, testMergeSubsets,
                            testMergeCircular, testAllToSingletonsAfterMergesRestoresDistinctSets,
                            testToPartitionUsesTemplateIndexType);

using GenericUnionFindTestTypes =
    ::testing::Types<index, std::uint32_t, std::uint16_t, std::int64_t, std::int32_t, std::int16_t>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestGenericUnionFind, GenericUnionFindGTest,
                               GenericUnionFindTestTypes, );

static_assert(std::same_as<UnionFind, GenericUnionFind<index>>);

} // namespace
} /* namespace NetworKit */
