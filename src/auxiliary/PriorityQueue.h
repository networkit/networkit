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
	std::map<Val, uint64_t> heapIndex;

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
	 * Removes key-value pair given by @a elem.
	 */
	virtual void remove(const ElemType& elem);

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

	/**
	 * Removes all elements from the PQ.
	 */
	virtual void clear();

	virtual void print() const;
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

	for (uint64_t i = 0; i < n; ++i) {
		heapIndex.insert(std::make_pair(elems[i].second, i));
	}

	pq = elems;

	if (n >= 2) {
		for (int64_t i = (n-2) / 2; i >= 0; --i) {
			TRACE("start heapifyDown for " << i);
			heapifyDown(i);
		}
	}

//	for (uint64_t i = 0; i < pq.size(); ++i) {
//		TRACE("binary heap elem  " << i << ": (" << pq[i].first << "," << pq[i].second << ")");
//	}
}


template<class Key, class Val>
Aux::PriorityQueue<Key, Val>::~PriorityQueue() {

}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::heapifyDown(uint64_t pos) {
	TRACE("PQ: heapifyDown");

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
		TRACE("exchange " << pos << " and " << smallest);
	      exchange(pos, smallest);
	      heapifyDown(smallest); // index smallest contains value of pos now
	}
}

template<class Key, class Val>
uint64_t Aux::PriorityQueue<Key, Val>::heapifyUp(uint64_t pos) {
	TRACE("PQ: heapifyUp");

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
	TRACE("PQ: insert");
	uint64_t pqPos = pq.size();
	pq.push_back(elem);

	// send element up as far as necessary to restore heap order
	TRACE("call heapifyUp for new elem");
	heapIndex.insert(std::make_pair(elem.second, pqPos));
	if (pq.size() > 1) {
		uint64_t newElemIdx = heapifyUp(pqPos);
		TRACE("heapifyUp reports new index " << newElemIdx);
		heapIndex[elem.second] = newElemIdx;
	}
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
	TRACE("PQ: decrease key");

	// find element
	uint64_t index = find(elem);
	bool found = (index != none);

	if (found) {
		// change key
		assert(elem.first < pq[index].first && elem.second == pq[index].second);
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
	TRACE("PQ: find");

	Val val = elem.second;
	if (heapIndex.count(val) > 0) {
		TRACE("PQ: find result: " << heapIndex.find(val)->second);
		return (heapIndex.find(val)->second); // TODO: avoid double access to map
	}
	else {
		TRACE("PQ: find result: " << none);
		return none;
	}
}

template<class Key, class Val>
inline void Aux::PriorityQueue<Key, Val>::exchange(uint64_t pos1, uint64_t pos2) {
	DEBUG("PQ: exchange, pos1: " << pos1 << ", pos2: " << pos2);

	assert(pos1 != none && pos2 != none);
	TRACE("heapIndex of size " << heapIndex.size() << ", swap indices " << pos1 << " and " << pos2);
	std::swap(pq[pos1], pq[pos2]);
//	TRACE("elements: " << pq[pos1].second << " and " << pq[pos2].second);
	heapIndex[pq[pos1].second] = pos1;
	heapIndex[pq[pos2].second] = pos2;
}

template<class Key, class Val>
inline uint64_t Aux::PriorityQueue<Key, Val>::size() const {
	return pq.size();
}

template<class Key, class Val>
void Aux::PriorityQueue<Key, Val>::remove(const ElemType& elem) {
	DEBUG("PQ: remove");

	// find element
	uint64_t pqSize = pq.size();
	bool found = (pqSize > 0);
	uint64_t elemIndex = find(elem);
	found &= (elemIndex != none);

	if (found) { // remove
		bool performExchange = (pqSize >= 2);

		if (performExchange) {
			// swap with last element, delete after swap
			uint64_t lastIndex = pqSize - 1;
			TRACE("PQ remove, pqSize: " << pqSize << ", elemIndex: " << elemIndex << ", lastIndex: " << lastIndex);
			this->exchange(elemIndex, lastIndex);
		}

		pq.pop_back();

		if (performExchange) {
			// put formerly last element into correct heap order
			this->heapifyDown(elemIndex);
		}

		heapIndex.erase(elem.second);
	}
//	else {
//		DEBUG("trying to remove non-existing element from PQ!");
//	}
}

template<class Key, class Val>
inline void Aux::PriorityQueue<Key, Val>::clear() {
	pq.clear();
	heapIndex.clear();
}

template<class Key, class Val>
inline void Aux::PriorityQueue<Key, Val>::print() const {
	uint64_t behindLineIndex = 1;
	uint64_t index = 0;
	uint64_t size = pq.size();

	if (size > 0) {
		while (behindLineIndex <= size) {
			for (uint64_t i = index; i < behindLineIndex; ++i) {
				std::cout << pq[i].first << "  ";
			}
			index = behindLineIndex;
			behindLineIndex = 2 * behindLineIndex + 1;
			std::cout << std::endl;
		}
	}
}

#endif /* PRIORITYQUEUE_H_ */
