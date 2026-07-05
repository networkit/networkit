/*
 * BucketPQ.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Henning
 */

#ifndef NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_
#define NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_

#include <cassert>
#include <concepts>
#include <cstdint>
#include <limits>
#include <list>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>

namespace Aux {

using index = NetworKit::index;
using count = NetworKit::count;

constexpr int64_t none = std::numeric_limits<int64_t>::max();

template <typename T>
concept SignedIntegral = std::integral<T> && std::is_signed_v<T>;

template <typename T>
concept IntegralValue = std::unsigned_integral<T> && !std::same_as<std::remove_cvref_t<T>, bool>;

/**
 * Addressable priority queue for values in the range [0,n) and
 * integer keys (= priorities) in the range [minPrio, maxPrio].
 * minPrio and maxPrio can be positive or negative, respectively with
 * the obvious constraint minPrio <= maxPrio.
 * Amortized constant running time for each operation.
 */
template <SignedIntegral KeyType = int64_t, IntegralValue ValueType = index>
class BucketPQ : public PrioQueue<KeyType, ValueType> {
private:
    using Bucket = std::list<ValueType>;
    using BucketIndex = index;

    static constexpr KeyType noneKey = std::numeric_limits<KeyType>::max();
    static constexpr ValueType noneValue = std::numeric_limits<ValueType>::max();
    static constexpr BucketIndex noneBucket = std::numeric_limits<BucketIndex>::max();

    std::vector<Bucket> buckets; // the actual buckets
    inline static Bucket dummyBucket{};
    inline static const Bucket::iterator invalidPtr = dummyBucket.end();

    struct OptionalIterator {
        bool valid;
        Bucket::iterator iter;

        void reset() {
            valid = false;
            iter = invalidPtr;
        }
        OptionalIterator() { reset(); }
        OptionalIterator(bool valid, Bucket::iterator iter) : valid(valid), iter(iter) {}
    };

    std::vector<OptionalIterator> nodePtr; // keeps track of node positions
    std::vector<BucketIndex> myBucket;           // keeps track of current bucket for each value
    KeyType currentMinKey;                 // current min key
    KeyType currentMaxKey;                 // current max key
    KeyType minAdmissibleKey;              // minimum admissible key
    KeyType maxAdmissibleKey;              // maximum admissible key
    count numElems;                        // number of elements stored
    KeyType offset;                        // offset from minAdmissibleKeys to 0

    /**
     * Constructor. Not to be used, only here for overriding.
     */
    BucketPQ(std::span<const KeyType>) {}

    /**
     * Constructor. Not to be used, only here for overriding.
     */
    BucketPQ(uint64_t) {}

    /**
     * Called from various constructors for initializing members.
     */
    void construct(uint64_t capacity) {
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

public:
    /**
     * Builds priority queue from the vector @a keys, values are indices
     * of @a keys.
     * @param[in] keys Vector of keys
     * @param[in] minAdmissibleKey Minimum admissible key
     * @param[in] maxAdmissibleKey Maximum admissible key
     */
    BucketPQ(std::span<const KeyType> keys, KeyType minAdmissibleKey, KeyType maxAdmissibleKey)
        : minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey) {
        construct(keys.size());

        // insert key-value pairs
        for (index i = 0; i < keys.size(); ++i) {
            if (keys[i] != noneKey) {
                insert(keys[i], static_cast<ValueType>(i));
            }
        }
    }

    /**
     * Builds priority queue of the specified capacity @a capacity.
     */
    BucketPQ(uint64_t capacity, KeyType minAdmissibleKey, KeyType maxAdmissibleKey)
        : minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey) {
        construct(capacity);
    }

    /**
     * Default destructor
     */
    ~BucketPQ() override = default;

    /**
     * Inserts key-value pair (@key, @value).
     */
    void insert(KeyType key, ValueType value) override {
        assert(minAdmissibleKey <= key && key <= maxAdmissibleKey);

        const auto valueIdx = static_cast<index>(value);
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

    /**
     * Returns the element on top of the PrioQ.
     */
    std::pair<KeyType, ValueType> getMin() {
        if (empty())
            return {noneKey, noneValue};
        else
            return {currentMinKey, buckets[currentMinKey + offset].front()};
    }

    /**
     * Removes the element with minimum key and returns the key-value pair.
     */
    std::pair<KeyType, ValueType> extractMin() override {
        if (empty())
            return {noneKey, noneValue};

        ValueType result = buckets[currentMinKey + offset].front();

        // store currentMinKey because remove(result) will change it
        KeyType oldMinKey = currentMinKey;
        remove(result);
        return {oldMinKey, result};
    }

    /**
     * Modifies entry with value @a value.
     * The entry is then set to @a newKey with the same value.
     * If the corresponding key is not present, the element will be inserted.
     */
    void changeKey(KeyType newKey, ValueType value) override {
        remove(value);
        insert(newKey, value);
    }

    /**
     * @return Number of elements in PQ.
     */
    uint64_t size() const override { return numElems; }

    /**
     * @return Whether or not the PQ is empty.
     */
    bool empty() const noexcept override { return numElems == 0; }

    /**
     * @return Whether or not the PQ contains the given element.
     */
    bool contains(const ValueType &value) const override {
        const auto valueIdx = static_cast<index>(value);
        return valueIdx < nodePtr.size() && nodePtr[valueIdx].valid;
    }

    /**
     * Removes key-value pair given by value @a val.
     */
    void remove(const ValueType &value) override {
        const auto valueIdx = static_cast<index>(value);
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

    /**
     * @return key to given value @val.
     */
    virtual KeyType getKey(const ValueType &val) {
        const auto valueIdx = static_cast<index>(val);
        return static_cast<KeyType>(myBucket[valueIdx]) - offset;
    }
};

} /* namespace Aux */
#endif // NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_
