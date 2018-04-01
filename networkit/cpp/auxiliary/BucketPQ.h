/*
 * BucketPQ.h
 *
 *  Created on: 02.03.2017
 *      Author: Henning
 */

#ifndef BUCKETPQ_H_
#define BUCKETPQ_H_

#include "Log.h"
#include "../Globals.h"
#include "PrioQueue.h"
#include <list>
#include <limits>

namespace Aux {

typedef NetworKit::index index;
typedef NetworKit::count count;
typedef std::list<index> Bucket;
constexpr int64_t none = std::numeric_limits<int64_t>::max();

/**
 * Addressable priority queue for values in the range [0,n) and
 * integer keys (= priorities) in the range [minPrio, maxPrio].
 * minPrio and maxPrio can be positive or negative, respectively with
 * the obvious constraint minPrio <= maxPrio.
 * Amortized constant running time for each operation.
 */
class BucketPQ: public PrioQueue<int64_t, index> {
private:
	std::vector<Bucket> buckets;			// the actual buckets
	std::vector<Bucket::iterator> nodePtr;	// keeps track of node positions
	std::vector<index> myBucket;			// keeps track of current bucket for each value
	int64_t currentMinKey;					// current min key
	int64_t currentMaxKey;					// current max key
	int64_t minAdmissibleKey;				// minimum admissible key
	int64_t maxAdmissibleKey;				// maximum admissible key
	count numElems;							// number of elements stored
	int64_t offset;							// offset from minAdmissibleKeys to 0

private:
	/**
	 * Constructor. Not to be used, only here for overriding.
	 */
	BucketPQ(const std::vector<int64_t>& keys) {}

	/**
	 * Constructor. Not to be used, only here for overriding.
	 */
	BucketPQ(uint64_t capacity) {}

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
	BucketPQ(const std::vector<int64_t>& keys, int64_t minAdmissibleKey, int64_t maxAdmissibleKey);

	/**
	* Builds priority queue of the specified capacity @a capacity.
	*/
	BucketPQ(uint64_t capacity, int64_t minAdmissibleKey, int64_t maxAdmissibleKey);

	/**
	 * Default destructor
	 */
	virtual ~BucketPQ() = default;

	/**
	 * Inserts key-value pair (@key, @value).
	 */
	virtual void insert(int64_t key, index value) override;

	/**
	 * Removes the element with minimum key and returns the key-value pair.
	 */
	virtual std::pair<int64_t, index> extractMin() override;

	/**
	 * Modifies entry with value @a value.
	 * The entry is then set to @a newKey with the same value.
	 * If the corresponding key is not present, the element will be inserted.
	 */
	virtual void changeKey(int64_t newKey, index value) override;

	/**
	 * @return Number of elements in PQ.
	 */
	virtual uint64_t size() const override;

	/**
	 * Removes key-value pair given by value @a val.
	 */
	virtual void remove(const index& val) override;

	/**
	 * @return key to given value @val.
	 */
	virtual int getKey(const index& val);
       	
};

} /* namespace Aux */
#endif /* BUCKETPQ_H_ */
