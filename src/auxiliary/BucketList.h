/*
 * BucketList.h
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#ifndef BUCKETLIST_H_
#define BUCKETLIST_H_

#include <vector>
#include <list>

namespace aux {

class BucketList {
private:
	std::vector<std::list<int64_t> > buckets;
	uint64_t firstList;

public:
	BucketList();
	virtual ~BucketList();

	/**
	 * Insert node entry @a v into list @a value.
	 */
	void insert(node v, count listNum);

	/**
	 * Access first value, move it to next smaller list or
	 * remove if in list 0.
	 */
	void decreaseFirstValue();
};

} /* namespace aux */
#endif /* BUCKETLIST_H_ */
