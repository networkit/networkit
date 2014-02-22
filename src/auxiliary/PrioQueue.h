/*
 * PrioQueue.h
 *
 *  Created on: 21.02.2014
 *      Author: Henning
 */

#ifndef PRIOQUEUE_H_
#define PRIOQUEUE_H_

#include <cassert>
#include <set>
#include <vector>
#include <limits>

namespace Aux {


/**
 * Priority queue with extract-min and decrease-key.
 * The type Val takes on integer values between 0 and n-1.
 * O(n log n) for construction, O(log n) for typical operations.
 */
template<class Key, class Val>
class PrioQueue {
	typedef std::pair<Key, Val> ElemType;

protected:
	std::set<ElemType> pqset;
	std::vector<Key> mapValToKey;

	/**
	 * Removes key-value pair given by @a elem.
	 */
	virtual void remove(const ElemType& elem);

public:
	/**
	 * Builds priority queue from the vector @a elems.
	 */
	PrioQueue(const std::vector<ElemType>& elems);

	/**
	 * Builds priority queue from the vector @a keys, values are indices
	 * in @a keys.
	 */
	PrioQueue(std::vector<Key>& keys);

	virtual ~PrioQueue() {}

	/**
	 * Inserts key-value pair stored in @a elem.
	 */
	virtual void insert(Key key, Val value);

	/**
	 * Removes the element with minimum key and returns it.
	 */
	virtual ElemType extractMin();

	/**
	 * Modifies entry with key @a elem.first and sets its value
	 * to @a elem.second. If the corresponding key is not present,
	 * the element will be inserted.
	 */
	virtual void decreaseKey(Key key, Val newValue);

	/**
	 * @return Number of elements in PQ.
	 */
	virtual uint64_t size() const;

	/**
	 * Removes all elements from the PQ.
	 */
	virtual void clear();
};

} /* namespace Aux */


template<class Key, class Val>
Aux::PrioQueue<Key, Val>::PrioQueue(const std::vector<ElemType>& elems) {
	mapValToKey.resize(elems.size());
	for (auto elem: elems) {
		insert(elem.first, elem.second);
	}
}

template<class Key, class Val>
Aux::PrioQueue<Key, Val>::PrioQueue(std::vector<Key>& keys) {
	mapValToKey.resize(keys.size());
	uint64_t index = 0;
	for (auto key: keys) {
		insert(key, index);
		++index;
	}
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::insert(Key key, Val value) {
	pqset.insert(std::make_pair(key, value));
	mapValToKey.at(value) = key;
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::remove(const ElemType& elem) {
	pqset.erase(elem);
	mapValToKey.at(elem.second) = none;
}

template<class Key, class Val>
std::pair<Key, Val> Aux::PrioQueue<Key, Val>::extractMin() {
	assert(pqset.size() > 0);
	ElemType elem = (* pqset.begin());
	remove(elem);
	mapValToKey.at(elem.second) = none;
	return elem;
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::decreaseKey(Key key, Val newValue) {
	// find and remove element with given key
	remove(std::make_pair(key, newValue));

	// insert element with new value
	insert(key, newValue);
}

template<class Key, class Val>
inline uint64_t Aux::PrioQueue<Key, Val>::size() const {
	return pqset.size();
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::clear() {
	pqset.clear();
	mapValToKey.clear();
}


#endif /* PRIOQUEUE_H_ */
