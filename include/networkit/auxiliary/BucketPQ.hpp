/*
 * BucketPQ.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Henning
 */

#ifndef NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_
#define NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_

#include <cassert>
#include <cstdint>
#include <limits>
#include <list>
#include <span>
#include <stdexcept>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>

namespace Aux {

using index = NetworKit::index;
using count = NetworKit::count;
using Bucket = std::list<index>;
constexpr int64_t none = std::numeric_limits<int64_t>::max();

/**
 * Addressable priority queue for values in the range [0,n) and
 * integer keys (= priorities) in the range [minPrio, maxPrio].
 * minPrio and maxPrio can be positive or negative, respectively with
 * the obvious constraint minPrio <= maxPrio.
 * Amortized constant running time for each operation.
 */
class BucketPQ : public PrioQueue<int64_t, index> {
private:
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
    std::vector<index> myBucket;           // keeps track of current bucket for each value
    int64_t currentMinKey;                 // current min key
    int64_t currentMaxKey;                 // current max key
    int64_t minAdmissibleKey;              // minimum admissible key
    int64_t maxAdmissibleKey;              // maximum admissible key
    count numElems;                        // number of elements stored
    int64_t offset;                        // offset from minAdmissibleKeys to 0

    /**
     * Constructor. Not to be used, only here for overriding.
     */
    BucketPQ(std::span<const int64_t>) {}

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
        myBucket.resize(capacity);
        currentMinKey = std::numeric_limits<int64_t>::max();
        currentMaxKey = std::numeric_limits<int64_t>::min();
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
    BucketPQ(std::span<const int64_t> keys, int64_t minAdmissibleKey, int64_t maxAdmissibleKey)
        : minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey) {
        construct(keys.size());

        // insert key-value pairs
        for (index i = 0; i < keys.size(); ++i) {
            if (keys[i] != none) {
                insert(keys[i], i);
            }
        }
    }

    /**
     * Builds priority queue of the specified capacity @a capacity.
     */
    BucketPQ(uint64_t capacity, int64_t minAdmissibleKey, int64_t maxAdmissibleKey)
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
    void insert(int64_t key, index value) override {
        assert(minAdmissibleKey <= key && key <= maxAdmissibleKey);
        assert(value < nodePtr.size());

        buckets[key + offset].push_front(value);
        nodePtr[value] = OptionalIterator{true, buckets[key + offset].begin()};
        myBucket[value] = key + offset;
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
    std::pair<int64_t, index> getMin() {
        if (empty())
            return {none, none};
        else
            return {currentMinKey, buckets[currentMinKey + offset].front()};
    }

    /**
     * Removes the element with minimum key and returns the key-value pair.
     */
    std::pair<int64_t, index> extractMin() override {
        if (empty())
            return {none, none};

        index result = buckets[currentMinKey + offset].front();

        // store currentMinKey because remove(result) will change it
        int64_t oldMinKey = currentMinKey;
        remove(result);
        return {oldMinKey, result};
    }

    /**
     * Modifies entry with value @a value.
     * The entry is then set to @a newKey with the same value.
     * If the corresponding key is not present, the element will be inserted.
     */
    void changeKey(int64_t newKey, index value) override {
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
    bool contains(const index &value) const override {
        return value < nodePtr.size() && nodePtr[value].valid;
    }

    /**
     * Removes key-value pair given by value @a val.
     */
    void remove(const index &value) override {
        assert(value < nodePtr.size());

        if (myBucket[value] != none) {
            // remove from appropriate bucket
            index key = myBucket[value];
            buckets[key].erase(nodePtr[value].iter);
            nodePtr[value].reset();
            myBucket[value] = none;
            --numElems;

            if (empty()) {
                // empty pq: reinit the current min/max pointers
                currentMinKey = std::numeric_limits<int64_t>::max();
                currentMaxKey = std::numeric_limits<int64_t>::min();
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
    virtual int64_t getKey(const index &val) {
        return static_cast<int64_t>(myBucket[val]) - offset;
    }
};

} /* namespace Aux */
#endif // NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_
