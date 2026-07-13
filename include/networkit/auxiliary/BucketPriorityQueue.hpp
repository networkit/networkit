/*
 * BucketPriorityQueue.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Henning
 */

#ifndef NETWORKIT_AUXILIARY_BUCKET_PRIORITY_QUEUE_HPP_
#define NETWORKIT_AUXILIARY_BUCKET_PRIORITY_QUEUE_HPP_

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>

namespace Aux {

using index = NetworKit::index;
using count = NetworKit::count;

template <typename T>
concept SignedIntegral = std::integral<T> && std::is_signed_v<T>;

template <typename T>
concept IntegralValue = std::integral<T> && !std::same_as<std::remove_cvref_t<T>, bool>;

/**
 * Addressable priority queue for values in the range [0,n) and
 * integer keys (= priorities) in the range [minPrio, maxPrio].
 * minPrio and maxPrio can be positive or negative, respectively with
 * the obvious constraint minPrio <= maxPrio.
 * Amortized constant running time for each operation.
 */
template <SignedIntegral KeyType = int64_t, IntegralValue ValueType = index>
class BucketPriorityQueue : public PrioQueue<KeyType, ValueType> {
private:
    using BucketIndex = index;

    static constexpr KeyType noneKey = std::numeric_limits<KeyType>::max();
    static constexpr ValueType noneValue = std::numeric_limits<ValueType>::max();
    static constexpr BucketIndex noneBucket = std::numeric_limits<BucketIndex>::max();

    std::vector<ValueType> bucketHead;
    std::vector<ValueType> previous;
    std::vector<ValueType> next;
    std::vector<BucketIndex> myBucket; // keeps track of current bucket for each value
    KeyType currentMinKey;             // current min key
    KeyType currentMaxKey;             // current max key
    KeyType minAdmissibleKey;          // minimum admissible key
    KeyType maxAdmissibleKey;          // maximum admissible key
    count numElems;                    // number of elements stored
    KeyType offset;                    // offset from minAdmissibleKeys to 0

    static index valueToIndex(ValueType value) {
        if constexpr (std::is_signed_v<ValueType>) {
            assert(value >= 0);
        }
        return static_cast<index>(value);
    }

    static void assureCapacityFitsValueType(uint64_t capacity) {
        if (capacity > static_cast<uint64_t>(noneValue)) {
            throw std::invalid_argument("capacity cannot be represented by ValueType");
        }
    }

    /**
     * Constructor. Not to be used, only here for overriding.
     */
    BucketPriorityQueue(std::span<const KeyType>) {}

    /**
     * Constructor. Not to be used, only here for overriding.
     */
    BucketPriorityQueue(uint64_t) {}

    /**
     * Called from various constructors for initializing members.
     */
    void construct(uint64_t capacity);

public:
    /**
     * Builds priority queue from the vector @a keys, values are indices
     * of @a keys.
     * @param[in] keys Vector of keys
     * @param[in] minAdmissibleKey Minimum admissible key
     * @param[in] maxAdmissibleKey Maximum admissible key
     */
    BucketPriorityQueue(std::span<const KeyType> keys, KeyType minAdmissibleKey,
                        KeyType maxAdmissibleKey);

    /**
     * Builds priority queue of the specified capacity @a capacity.
     */
    BucketPriorityQueue(uint64_t capacity, KeyType minAdmissibleKey, KeyType maxAdmissibleKey);

    /**
     * Default destructor
     */
    ~BucketPriorityQueue() override = default;

    /**
     * Inserts key-value pair (@key, @value).
     */
    void insert(KeyType key, ValueType value) override;

    /**
     * Returns the element on top of the PrioQ.
     */
    std::pair<KeyType, ValueType> getMin();

    /**
     * Removes the element with minimum key and returns the key-value pair.
     */
    std::pair<KeyType, ValueType> extractMin() override;

    /**
     * Modifies entry with value @a value.
     * The entry is then set to @a newKey with the same value.
     * If the corresponding key is not present, the element will be inserted.
     */
    void changeKey(KeyType newKey, ValueType value) override;

    /**
     * @return Number of elements in PQ.
     */
    uint64_t size() const override;

    /**
     * @return Whether or not the PQ is empty.
     */
    bool empty() const noexcept override;

    /**
     * @return Whether or not the PQ contains the given element.
     */
    bool contains(const ValueType &value) const override;

    /**
     * Removes key-value pair given by value @a val.
     */
    void remove(const ValueType &value) override;

    /**
     * Removes all elements from the priority queue.
     */
    void clear() override;

    /**
     * @return key to given value @val.
     */
    virtual KeyType getKey(const ValueType &val);
};

using BucketPQ = BucketPriorityQueue<int64_t, index>;

} /* namespace Aux */

#include <networkit/auxiliary/BucketPriorityQueueImpl.hpp>

#endif // NETWORKIT_AUXILIARY_BUCKET_PRIORITY_QUEUE_HPP_
