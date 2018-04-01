/*
 * BucketPQ.cpp
 *
 *  Created on: 03.03.2017
 *      Author: Henning
 */

#include "BucketPQ.h"

namespace Aux {

BucketPQ::BucketPQ(const std::vector<int64_t>& keys, int64_t minAdmissibleKey, int64_t maxAdmissibleKey):
		minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey)
{
	construct(keys.size());

	// insert key-value pairs
	for (index i = 0; i < keys.size(); ++i) {
		if (keys[i] != none) {
			insert(keys[i], i);
		}
	}
}

BucketPQ::BucketPQ(uint64_t capacity, int64_t minAdmissibleKey, int64_t maxAdmissibleKey):
		minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey)
{
	construct(capacity);
}

void BucketPQ::construct(uint64_t capacity) {
	// check range of keys
	if (minAdmissibleKey > maxAdmissibleKey) {
		throw std::invalid_argument("minAdmissibleKey cannot be larger than maxAdmissibleKey");
	}

	// init
	buckets.resize(maxAdmissibleKey-minAdmissibleKey+1);
	nodePtr.resize(capacity);
	myBucket.resize(capacity);
	currentMinKey = std::numeric_limits<int64_t>::max();
	currentMaxKey = std::numeric_limits<int64_t>::min();
	numElems = 0;

	offset = -minAdmissibleKey;
}

void BucketPQ::insert(int64_t key, index value) {
	assert(minAdmissibleKey <= key && key <= maxAdmissibleKey);
	assert(0 <= value && value < nodePtr.size());

	buckets[key+offset].push_front(value);
	nodePtr[value] = buckets[key+offset].begin();
	myBucket[value] = key+offset;
	++numElems;

	// bookkeeping
	if (key < currentMinKey) {
		currentMinKey = key;
	}
	if (key > currentMaxKey) {
		currentMaxKey = key;
	}
}

void BucketPQ::remove(const index& value) {
	assert(0 <= value && value < nodePtr.size());

	if (myBucket[value] != none) {
		// remove from appropriate bucket
		index key = myBucket[value];
		buckets[key].erase(nodePtr[value]);
		myBucket[value] = none;
		--numElems;

		if (size() == 0) {
			// empty pq: reinit the current min/max pointers
			currentMinKey = std::numeric_limits<int64_t>::max();
			currentMaxKey = std::numeric_limits<int64_t>::min();
		}
		else {
			// adjust max pointer if necessary
			while (buckets[currentMaxKey+offset].empty() && currentMaxKey > currentMinKey) {
				--currentMaxKey;
			}

			// adjust min pointer if necessary
			while (buckets[currentMinKey+offset].empty() && currentMinKey < currentMaxKey) {
				++currentMinKey;
			}
		}
	}
}

std::pair<int64_t, index> BucketPQ::extractMin() {
	if (size() == 0) {
		return std::make_pair(none, NetworKit::none);
	}
	else {
		assert(! buckets[currentMinKey+offset].empty());
		index result = buckets[currentMinKey+offset].front();

		// store currentMinKey because remove(result) will change it
		int64_t oldMinKey = currentMinKey;
		remove(result);
		return std::make_pair(oldMinKey, result);
	}
}

void BucketPQ::changeKey(int64_t newKey, index value) {
	remove(value);
	insert(newKey, value);
}

uint64_t BucketPQ::size() const {
	return numElems;
}

int BucketPQ::getKey(const index& val) {
	int key =myBucket[val]-offset;
	return key;
}
} // namespace
