/*
 * BucketPriorityQueueImpl.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Henning
 */

#ifndef NETWORKIT_AUXILIARY_BUCKET_PRIORITY_QUEUE_IMPL_HPP_
#define NETWORKIT_AUXILIARY_BUCKET_PRIORITY_QUEUE_IMPL_HPP_

#include <networkit/auxiliary/BucketPriorityQueue.hpp>

namespace Aux {

template <SignedIntegral KeyType, IntegralValue ValueType>
void BucketPriorityQueue<KeyType, ValueType>::construct(uint64_t capacity) {
    // check range of keys
    if (minAdmissibleKey > maxAdmissibleKey) {
        throw std::invalid_argument("minAdmissibleKey cannot be larger than maxAdmissibleKey");
    }

    // init
    buckets.resize(maxAdmissibleKey - minAdmissibleKey + 1);
    nodePtr.resize(capacity);
    myBucket.assign(capacity, noneBucket);
    currentMinKey = noneKey;
    currentMaxKey = std::numeric_limits<KeyType>::min();
    numElems = 0;

    offset = -minAdmissibleKey;
}

template <SignedIntegral KeyType, IntegralValue ValueType>
BucketPriorityQueue<KeyType, ValueType>::BucketPriorityQueue(std::span<const KeyType> keys,
                                                             KeyType minAdmissibleKey,
                                                             KeyType maxAdmissibleKey)
    : minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey) {
    construct(keys.size());

    // insert key-value pairs
    for (index i = 0; i < keys.size(); ++i) {
        if (keys[i] != noneKey) {
            insert(keys[i], static_cast<ValueType>(i));
        }
    }
}

template <SignedIntegral KeyType, IntegralValue ValueType>
BucketPriorityQueue<KeyType, ValueType>::BucketPriorityQueue(uint64_t capacity,
                                                             KeyType minAdmissibleKey,
                                                             KeyType maxAdmissibleKey)
    : minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey) {
    construct(capacity);
}

template <SignedIntegral KeyType, IntegralValue ValueType>
void BucketPriorityQueue<KeyType, ValueType>::insert(KeyType key, ValueType value) {
    assert(minAdmissibleKey <= key && key <= maxAdmissibleKey);

    const auto valueIdx = static_cast<index>(value);
    if constexpr (std::is_signed_v<ValueType>) {
        assert(value >= 0);
    }
    assert(valueIdx < nodePtr.size());

    const auto bucketIdx = static_cast<BucketIndex>(key + offset);
    buckets[bucketIdx].push_front(value);
    nodePtr[valueIdx] = OptionalIterator{true, buckets[bucketIdx].begin()};
    myBucket[valueIdx] = bucketIdx;
    ++numElems;

    // bookkeeping
    if (key < currentMinKey) {
        currentMinKey = key;
    }
    if (key > currentMaxKey) {
        currentMaxKey = key;
    }
}

template <SignedIntegral KeyType, IntegralValue ValueType>
std::pair<KeyType, ValueType> BucketPriorityQueue<KeyType, ValueType>::getMin() {
    if (empty())
        return {noneKey, noneValue};
    else
        return {currentMinKey, buckets[currentMinKey + offset].front()};
}

template <SignedIntegral KeyType, IntegralValue ValueType>
std::pair<KeyType, ValueType> BucketPriorityQueue<KeyType, ValueType>::extractMin() {
    if (empty())
        return {noneKey, noneValue};

    ValueType result = buckets[currentMinKey + offset].front();

    // store currentMinKey because remove(result) will change it
    KeyType oldMinKey = currentMinKey;
    remove(result);
    return {oldMinKey, result};
}

template <SignedIntegral KeyType, IntegralValue ValueType>
void BucketPriorityQueue<KeyType, ValueType>::changeKey(KeyType newKey, ValueType value) {
    remove(value);
    insert(newKey, value);
}

template <SignedIntegral KeyType, IntegralValue ValueType>
uint64_t BucketPriorityQueue<KeyType, ValueType>::size() const {
    return numElems;
}

template <SignedIntegral KeyType, IntegralValue ValueType>
bool BucketPriorityQueue<KeyType, ValueType>::empty() const noexcept {
    return numElems == 0;
}

template <SignedIntegral KeyType, IntegralValue ValueType>
bool BucketPriorityQueue<KeyType, ValueType>::contains(const ValueType &value) const {
    const auto valueIdx = static_cast<index>(value);
    if constexpr (std::is_signed_v<ValueType>) {
        assert(value >= 0);
    }
    return valueIdx < nodePtr.size() && nodePtr[valueIdx].valid;
}

template <SignedIntegral KeyType, IntegralValue ValueType>
void BucketPriorityQueue<KeyType, ValueType>::remove(const ValueType &value) {
    const auto valueIdx = static_cast<index>(value);
    if constexpr (std::is_signed_v<ValueType>) {
        assert(value >= 0);
    }
    assert(valueIdx < nodePtr.size());

    if (myBucket[valueIdx] != noneBucket) {
        // remove from appropriate bucket
        BucketIndex bucketIdx = myBucket[valueIdx];
        buckets[bucketIdx].erase(nodePtr[valueIdx].iter);
        nodePtr[valueIdx].reset();
        myBucket[valueIdx] = noneBucket;
        --numElems;

        if (empty()) {
            // empty pq: reinit the current min/max pointers
            currentMinKey = noneKey;
            currentMaxKey = std::numeric_limits<KeyType>::min();
        } else {
            // adjust max pointer if necessary
            while (buckets[currentMaxKey + offset].empty() && currentMaxKey > currentMinKey) {
                --currentMaxKey;
            }

            // adjust min pointer if necessary
            while (buckets[currentMinKey + offset].empty() && currentMinKey < currentMaxKey) {
                ++currentMinKey;
            }
        }
    }
}

template <SignedIntegral KeyType, IntegralValue ValueType>
KeyType BucketPriorityQueue<KeyType, ValueType>::getKey(const ValueType &val) {
    const auto valueIdx = static_cast<index>(val);
    if constexpr (std::is_signed_v<ValueType>) {
        assert(val >= 0);
    }
    return static_cast<KeyType>(myBucket[valueIdx]) - offset;
}

} /* namespace Aux */

#endif // NETWORKIT_AUXILIARY_BUCKET_PRIORITY_QUEUE_IMPL_HPP_
