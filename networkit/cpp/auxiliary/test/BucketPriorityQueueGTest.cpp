#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include <networkit/auxiliary/BucketPriorityQueue.hpp>

namespace NetworKit {
namespace {

using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::SizeIs;

template <typename KeyT, typename ValueT>
struct BucketPQConfig {
    using KeyType = KeyT;
    using ValueType = ValueT;
};

template <typename TestT>
class BucketPriorityQueueGTest : public ::testing::Test {
public:
    using KeyType = typename TestT::KeyType;
    using ValueType = typename TestT::ValueType;
};

TYPED_TEST_SUITE_P(BucketPriorityQueueGTest);

TYPED_TEST_P(BucketPriorityQueueGTest, testConstructFromKeysSkipsNone) {
    using KeyType = typename TestFixture::KeyType;
    using ValueType = typename TestFixture::ValueType;

    std::vector<KeyType> keys = {KeyType{3}, std::numeric_limits<KeyType>::max(), KeyType{-2},
                                 KeyType{1}};
    Aux::BucketPriorityQueue<KeyType, ValueType> prioQ(keys, KeyType{-5}, KeyType{5});

    EXPECT_THAT(prioQ, SizeIs(3));
    EXPECT_TRUE(prioQ.contains(ValueType{0}));
    EXPECT_FALSE(prioQ.contains(ValueType{1}));
    EXPECT_TRUE(prioQ.contains(ValueType{2}));
    EXPECT_TRUE(prioQ.contains(ValueType{3}));

    EXPECT_EQ(prioQ.getKey(ValueType{0}), KeyType{3});
    EXPECT_EQ(prioQ.getKey(ValueType{2}), KeyType{-2});
    EXPECT_EQ(prioQ.getKey(ValueType{3}), KeyType{1});

    EXPECT_EQ(prioQ.extractMin(), (std::pair<KeyType, ValueType>{KeyType{-2}, ValueType{2}}));
    EXPECT_EQ(prioQ.extractMin(), (std::pair<KeyType, ValueType>{KeyType{1}, ValueType{3}}));
    EXPECT_EQ(prioQ.extractMin(), (std::pair<KeyType, ValueType>{KeyType{3}, ValueType{0}}));
    EXPECT_THAT(prioQ, IsEmpty());
}

TYPED_TEST_P(BucketPriorityQueueGTest, testChangeKeyAndRemove) {
    using KeyType = typename TestFixture::KeyType;
    using ValueType = typename TestFixture::ValueType;

    Aux::BucketPriorityQueue<KeyType, ValueType> prioQ(5, KeyType{-10}, KeyType{10});

    prioQ.insert(KeyType{4}, ValueType{0});
    prioQ.insert(KeyType{-3}, ValueType{1});
    prioQ.insert(KeyType{2}, ValueType{2});

    EXPECT_TRUE(prioQ.contains(ValueType{0}));
    EXPECT_TRUE(prioQ.contains(ValueType{1}));
    EXPECT_TRUE(prioQ.contains(ValueType{2}));
    EXPECT_EQ(prioQ.getMin(), (std::pair<KeyType, ValueType>{KeyType{-3}, ValueType{1}}));

    prioQ.changeKey(KeyType{-5}, ValueType{2});
    EXPECT_EQ(prioQ.getKey(ValueType{2}), KeyType{-5});
    EXPECT_EQ(prioQ.extractMin(), (std::pair<KeyType, ValueType>{KeyType{-5}, ValueType{2}}));

    prioQ.remove(ValueType{1});
    EXPECT_FALSE(prioQ.contains(ValueType{1}));
    EXPECT_EQ(prioQ.extractMin(), (std::pair<KeyType, ValueType>{KeyType{4}, ValueType{0}}));
    EXPECT_THAT(prioQ, IsEmpty());
}

TYPED_TEST_P(BucketPriorityQueueGTest, testRemoveElementsWithSameKey) {
    using KeyType = typename TestFixture::KeyType;
    using ValueType = typename TestFixture::ValueType;

    Aux::BucketPriorityQueue<KeyType, ValueType> prioQ(6, KeyType{-2}, KeyType{2});

    prioQ.insert(KeyType{1}, ValueType{0});
    prioQ.insert(KeyType{1}, ValueType{1});
    prioQ.insert(KeyType{1}, ValueType{2});
    prioQ.insert(KeyType{-1}, ValueType{3});

    prioQ.remove(ValueType{1});
    EXPECT_FALSE(prioQ.contains(ValueType{1}));
    EXPECT_TRUE(prioQ.contains(ValueType{0}));
    EXPECT_TRUE(prioQ.contains(ValueType{2}));

    EXPECT_EQ(prioQ.extractMin(), (std::pair<KeyType, ValueType>{KeyType{-1}, ValueType{3}}));

    prioQ.remove(ValueType{2});
    EXPECT_FALSE(prioQ.contains(ValueType{2}));
    EXPECT_EQ(prioQ.extractMin(), (std::pair<KeyType, ValueType>{KeyType{1}, ValueType{0}}));
    EXPECT_THAT(prioQ, IsEmpty());
}

TYPED_TEST_P(BucketPriorityQueueGTest, testEmptySentinelsUseTemplateTypes) {
    using KeyType = typename TestFixture::KeyType;
    using ValueType = typename TestFixture::ValueType;

    Aux::BucketPriorityQueue<KeyType, ValueType> prioQ(4, KeyType{-2}, KeyType{2});

    const std::pair<KeyType, ValueType> emptySentinel{std::numeric_limits<KeyType>::max(),
                                                      std::numeric_limits<ValueType>::max()};
    EXPECT_THAT(prioQ, IsEmpty());
    EXPECT_EQ(prioQ.getMin(), emptySentinel);
    EXPECT_EQ(prioQ.extractMin(), emptySentinel);
}

TYPED_TEST_P(BucketPriorityQueueGTest, testClearResetsBucketState) {
    using KeyType = typename TestFixture::KeyType;
    using ValueType = typename TestFixture::ValueType;

    Aux::BucketPriorityQueue<KeyType, ValueType> prioQ(5, KeyType{-2}, KeyType{4});

    prioQ.insert(KeyType{0}, ValueType{0});
    prioQ.insert(KeyType{2}, ValueType{1});
    prioQ.insert(KeyType{-1}, ValueType{2});

    prioQ.clear();

    EXPECT_THAT(prioQ, IsEmpty());
    EXPECT_FALSE(prioQ.contains(ValueType{0}));
    EXPECT_THAT(prioQ.getMin(),
                Pair(std::numeric_limits<KeyType>::max(), std::numeric_limits<ValueType>::max()));

    prioQ.insert(KeyType{4}, ValueType{3});
    EXPECT_THAT(prioQ.extractMin(), Pair(KeyType{4}, ValueType{3}));
    EXPECT_THAT(prioQ, IsEmpty());
}

REGISTER_TYPED_TEST_SUITE_P(BucketPriorityQueueGTest, testConstructFromKeysSkipsNone,
                            testChangeKeyAndRemove, testEmptySentinelsUseTemplateTypes,
                            testClearResetsBucketState, testRemoveElementsWithSameKey);

using BucketPriorityQueueTestTypes =
    ::testing::Types<BucketPQConfig<int64_t, index>, BucketPQConfig<int32_t, uint32_t>,
                     BucketPQConfig<int16_t, uint16_t>, BucketPQConfig<int64_t, int>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestBucketPQ, BucketPriorityQueueGTest,
                               BucketPriorityQueueTestTypes, );

} // namespace
} // namespace NetworKit
