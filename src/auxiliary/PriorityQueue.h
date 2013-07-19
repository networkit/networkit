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

const uint64_t none = std::numeric_limits<uint64_t>::max();

/**
 * Priority queue with extract-min and decrease-key.
 * Realized as binary heap with an additional array pointing
 * to heap positions.
 * The values need to be distinct integers in the range [0..size-1],
 * where size is the number of elements for which build was called.
 */
template<class Key, class Val>
class PriorityQueue {
	typedef std::pair<Key, Val> ElemType;
	typedef std::vector<ElemType> Heap;

protected:
	Heap pq;
	std::vector<int64_t> heapIndex;

	/*** MIN-HEAP OPERATIONS ***/

	virtual void exchange(uint64_t pos1, uint64_t pos2);

	/**
	 * Restore heap property by moving element as far down as necessary.
	 */
	virtual void heapifyDown(uint64_t pos);

	/**
	 * Restore heap property by moving element as far up as necessary.
	 */
	virtual uint64_t heapifyUp(uint64_t pos);

	/**
	 * @return Left child of heap node at position @a pos.
	 */
	virtual uint64_t left(uint64_t pos) const;

	/**
	 * @return Right child of heap node at position @a pos.
	 */
	virtual uint64_t right(uint64_t pos) const;

	/**
	 * @return Parent of heap node at position @a pos.
	 */
	virtual uint64_t parent(uint64_t pos) const;

	/**
	 * @return Index of heap node with key @a elem.first.
	 */
	virtual uint64_t find(const ElemType& elem) const;

	/**
	 * Builds priority queue from the vector @a elems.
	 */
	virtual void init(const std::vector<ElemType>& elems);

public:
	/**
	 * Builds priority queue from the vector @a elems.
	 */
	PriorityQueue(const std::vector<ElemType>& elems);

	/**
	 * Builds priority queue from the vector @a keys, values are indices
	 * in @a keys.
	 */
	PriorityQueue(std::vector<Key>& keys);

	virtual ~PriorityQueue();

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

	/**
	 * @return Number of elements in PQ.
	 */
	virtual uint64_t size() const;
};

} /* namespace Aux */



template<class Key, class Val>
Aux::PriorityQueue<Key, Val>::PriorityQueue(const std::vector<ElemType>& elems) {
	init(elems);
}

template<class Key, class Val>
Aux::PriorityQueue<Key, Val>::PriorityQueue(std::vector<Key>& keys) {
	std::vector<ElemType> elems(keys.size());

	for (uint64_t i = 0; i < keys.size(); ++i) {
		elems[i] = std::make_pair(keys[i], i);
	}

	init(elems);
}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::init(const std::vector<ElemType>& elems) {
	uint64_t n = elems.size();

	heapIndex.resize(n);
	for (uint64_t i = 0; i < n; ++i) {
		heapIndex[i] = i;
	}

	pq = elems;
	for (int64_t i = (pq.size()-2) / 2; i >= 0; --i) {
		TRACE("start heapifyDown for " << i);
		heapifyDown(i);
	}

	for (uint64_t i = 0; i < pq.size(); ++i) {
		TRACE("binary heap elem  " << i << ": (" << pq[i].first << "," << pq[i].second << ")");
	}
}


template<class Key, class Val>
Aux::PriorityQueue<Key, Val>::~PriorityQueue() {

}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::heapifyDown(uint64_t pos) {
	// check if children are smaller
	uint64_t le = left(pos);
	uint64_t ri = right(pos);
	uint64_t smallest = pos;

	TRACE("pos: " << pos);
	TRACE("left: " << le);
	TRACE("right: " << ri);

	// test left child
	if ((le < pq.size()) && (pq[le].first < pq[pos].first)) {
		smallest = le;
	}

	// test right child
	if ((ri < pq.size()) && (pq[ri].first < pq[smallest].first)) {
		smallest = ri;
	}

	// move element down as long as necessary to restore heap order
	if (smallest != pos) {
	      exchange(pos, smallest);
	      heapifyDown(smallest); // index smallest contains value of pos now
	}
}

template<class Key, class Val>
uint64_t Aux::PriorityQueue<Key, Val>::heapifyUp(uint64_t pos) {
	// send element up as far as necessary to restore heap order
	uint64_t index = pos;
	uint64_t parentIdx = parent(index);
	while ((index > 0) && (pq[index].first < pq[parentIdx].first)) {
		this->exchange(index, parentIdx);
		index = parent(index);
		parentIdx = parent(index);
	}
	return index;
}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::insert(ElemType elem) {
	assert(elem.second >= pq.size());
	pq.push_back(elem);

	// send element up as far as necessary to restore heap order
	uint64_t newElemIdx = heapifyUp(pq.size() - 1);
	heapIndex.push_back(newElemIdx);
}

template<class Key, class Val>
std::pair<Key, Val> Aux::PriorityQueue<Key, Val>::extractMin() {
	assert(!pq.empty());
	ElemType elem = pq[0];

	// swap first and last element
	uint64_t n = pq.size();
	std::swap(pq[0], pq[n - 1]);

	// remove last element and restore heap property
	pq.pop_back();
	heapifyDown(0);

	return elem;
}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::decreaseKey(ElemType elem) {
	// find element
	uint64_t index = find(elem);
	bool found = index != none;

	if (found) {
		// change key
		ElemType& pqElem = pq[index];
		assert(elem.first < pqElem.first && elem.second == pqElem.second);
		pq[index].first = elem.first;

		// send element up as far as necessary to restore heap order
		this->heapifyUp(index);
	}
	else {	// if not found: insert elem
		WARN("decrease key on non-existing element with key " << elem.first
				<< ", insert instead!");
		this->insert(elem);
	}
}

template<class Key, class Val>
inline uint64_t Aux::PriorityQueue<Key, Val>::left(uint64_t pos) const {
	return 2*pos + 1;
}

template<class Key, class Val>
inline uint64_t Aux::PriorityQueue<Key, Val>::right(uint64_t pos) const {
	return 2*pos + 2;
}

template<class Key, class Val>
inline uint64_t Aux::PriorityQueue<Key, Val>::parent(uint64_t pos) const {
	return (pos - 1) / 2;
}

template<class Key, class Val>
inline uint64_t Aux::PriorityQueue<Key, Val>::find(const ElemType& elem) const {
	if (! (elem.second >= 0 && elem.second < heapIndex.size())) {
		WARN("find tries to access heap out of bounds");
		return none;
	}
	return heapIndex[elem.second];
}

template<class Key, class Val>
inline void Aux::PriorityQueue<Key, Val>::exchange(uint64_t pos1, uint64_t pos2) {
	std::swap(pq[pos1], pq[pos2]);
	TRACE("swap indices " << pq[pos1].second << " and " << pq[pos2].second);
	std::swap(heapIndex[pq[pos1].second], heapIndex[pq[pos2].second]);
}

template<class Key, class Val>
inline uint64_t Aux::PriorityQueue<Key, Val>::size() const {
	return pq.size();
}


#endif /* PRIORITYQUEUE_H_ */
