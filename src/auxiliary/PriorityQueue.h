/*
 * PriorityQueue.h
 *
 *  Created on: Jul 18, 2013
 *      Author: Henning
 */

#ifndef PRIORITYQUEUE_H_
#define PRIORITYQUEUE_H_

#include <algorithm>
#include <vector>

namespace Aux {

/**
 * Priority queue with extract-min and decrease-key.
 * Realized as binary heap with an additional array pointing
 * to heap positions.
 */
template<class Key, class Val>
class PriorityQueue {
	typedef std::pair<Key, Val> ElemType;
	typedef std::vector<ElemType> Heap;

protected:
	Heap pq;
	std::vector<int64_t> heapIndex;

	/*** MIN-HEAP OPERATIONS ***/

	void heapify();


public:
	PriorityQueue();
	virtual ~PriorityQueue();

	/**
	 * Builds priority queue from the vector @a elems.
	 */
	virtual void build(std::vector<ElemType>& elems);

	/**
	 * Inserts key-value pair stored in @a elem.
	 */
	virtual void insert(ElemType elem);

	/**
	 * Removes the element with minimum key and returns it.
	 */
	virtual ElemType extractMin();

	/**
	 * Modifies entry with key @a elem.first and sets its value
	 * to @a elem.second. If the corresponding key is not present,
	 * the element will be inserted.
	 */
	virtual void decreaseKey(ElemType elem);
};

} /* namespace Aux */


template<class Key, class Val>
Aux::PriorityQueue<Key, Val>::PriorityQueue() {

}

template<class Key, class Val>
Aux::PriorityQueue<Key, Val>::~PriorityQueue() {

}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::heapify() {
	// TODO
}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::build(std::vector<ElemType>& elems) {
	pq = elems;
	std::make_heap(pq.begin(), pq.end());

	uint64_t n = pq.size();
	heapIndex.resize(n);
	for (uint64_t i = 0; i < n; ++i) {
		uint64_t index = pq[i].second;
		heapIndex[index] = i;
	}
}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::insert(ElemType elem) {
	// TODO
}

template<class Key, class Val>
std::pair<Key, Val> Aux::PriorityQueue<Key, Val>::extractMin() {
	// TODO
}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::decreaseKey(ElemType elem) {
	// TODO
}

#endif /* PRIORITYQUEUE_H_ */
