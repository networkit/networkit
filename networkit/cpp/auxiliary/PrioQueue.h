/*
 * PrioQueue.h
 *
 *  Created on: 02.03.2017
 *      Author: Henning
 */

#ifndef PRIORITYQUEUE_H_
#define PRIORITYQUEUE_H_

#include <set>
#include <vector>

#include "../auxiliary/Log.h"

namespace Aux {

/**
 * Priority queue with extract-min and decrease-key.
 * The type Val takes on integer values between 0 and n-1.
 * O(n log n) for construction, O(log n) for typical operations.
 */
template<class Key, class Value>
class PrioQueue {
private:
	typedef std::pair<Key, Value> ElemType;

	std::set<ElemType> pqset; // TODO: would std::map work and simplify things?
	std::vector<Key> mapValToKey;

	const Key undefined = std::numeric_limits<Key>::max(); // TODO: make static


protected:

	/**
	 * Default constructor without functionality. Only here for derived classes!
	 */
	PrioQueue() = default;

	/**
	 * Removes key-value pair given by @a elem.
	 */
	virtual void remove(const ElemType& elem);

	/**
	 * @return current content of queue
	 */
	virtual std::set<std::pair<Key, Value>> content() const;


public:

	/**
	 * Builds priority queue from the vector @a keys, values are indices
	 * of @a keys.
	 */
	PrioQueue(const std::vector<Key>& keys);

	/**
	* Builds priority queue of the specified capacity @a capacity.
	*/
	PrioQueue(uint64_t capacity);

	/**
	 * Default destructor
	 */
	virtual ~PrioQueue() = default;

	/**
	 * Inserts key-value pair stored in @a elem.
	 */
	virtual void insert(Key key, Value value);

	/**
	 * Returns the @a n-th element in the priority queue.
	 */
	virtual ElemType peekMin(size_t n = 0);

	/**
	 * Removes the element with minimum key and returns it.
	 */
	virtual ElemType extractMin();

	/**
	 * Modifies entry with value @a value.
	 * The entry is then set to @a newKey with the same value.
	 * If the corresponding key is not present, the element will be inserted.
	 */
	virtual void changeKey(Key newKey, Value value);

	/**
	 * @return Number of elements in PQ.
	 */
	virtual uint64_t size() const;

	/**
	 * Removes key-value pair given by value @a val.
	 */
	virtual void remove(const Value& val);

	/**
	 * Removes all the elements from the priority queue.
	 */
	virtual void clear();

	/**
	 * Iterates over all the elements of the priority queue and call @a handle
	 * (lambda closure).
	 */
	template<typename L>
	void forElements(L handle) const;

	/**
	 * Iterates over all the elements of the priority queue while the condition
	 * is satisfied and call @a handle (lambda closure).
	 */
	template<typename C, typename L>
	void forElementsWhile(C condition, L handle) const;

	/**
	 * DEBUGGING
	 */
	virtual void print() {
		DEBUG("num entries: ", mapValToKey.size());
		for (uint64_t i = 0; i < mapValToKey.size(); ++i) {
			DEBUG("key: ", mapValToKey[i], ", val: ", i, "\n");
		}
	}
};


template<class Key, class Value>
Aux::PrioQueue<Key, Value>::PrioQueue(const std::vector<Key>& keys) {
	mapValToKey.resize(keys.size());
	uint64_t index = 0;
	for (auto key: keys) {
		insert(key, index);
		++index;
	}
}

template<class Key, class Value>
Aux::PrioQueue<Key, Value>::PrioQueue(uint64_t capacity) {
	mapValToKey.resize(capacity);
}

template<class Key, class Value>
inline void Aux::PrioQueue<Key, Value>::insert(Key key, Value value) {
	if (value >= mapValToKey.size()) {
		uint64_t doubledSize = 2 * mapValToKey.size();
		assert(value < doubledSize);
		mapValToKey.resize(doubledSize);
	}
	pqset.insert(std::make_pair(key, value));
	mapValToKey.at(value) = key;
}

template<class Key, class Value>
inline void Aux::PrioQueue<Key, Value>::remove(const ElemType& elem) {
	remove(elem.second);
}

template<class Key, class Value>
inline void Aux::PrioQueue<Key, Value>::remove(const Value& val) {
	Key key = mapValToKey.at(val);
//	DEBUG("key: ", key);
	pqset.erase(std::make_pair(key, val));
	mapValToKey.at(val) = undefined;
}

template<class Key, class Value>
std::pair<Key, Value> Aux::PrioQueue<Key, Value>::peekMin(size_t n) {
	assert(pqset.size() > n);
	ElemType elem = *std::next(pqset.begin(), n);
	return elem;
}

template<class Key, class Value>
std::pair<Key, Value> Aux::PrioQueue<Key, Value>::extractMin() {
	assert(pqset.size() > 0);
	ElemType elem = (* pqset.begin());
	remove(elem);
	return elem;
}

template<class Key, class Value>
inline void Aux::PrioQueue<Key, Value>::changeKey(Key newKey, Value value) {
	// find and remove element with given key
	remove(value);

	// insert element with new value
	insert(newKey, value);
}

template<class Key, class Value>
inline uint64_t Aux::PrioQueue<Key, Value>::size() const {
	return pqset.size();
}

template<class Key, class Value>
inline std::set<std::pair<Key, Value>> Aux::PrioQueue<Key, Value>::content() const {
	return pqset;
}

template<class Key, class Value>
inline void Aux::PrioQueue<Key, Value>::clear() {
	pqset.clear();
	auto capacity = mapValToKey.size();
	mapValToKey.clear();
	mapValToKey.resize(capacity);
}

template<class Key, class Value>
template<typename L>
inline void Aux::PrioQueue<Key, Value>::forElements(L handle) const {
	for (auto it = pqset.begin(); it != pqset.end(); ++it) {
		ElemType elem = *it;
		handle(elem.first, elem.second);
	}
}

template<class Key, class Value>
template<typename C, typename L>
inline void Aux::PrioQueue<Key, Value>::forElementsWhile(C condition,
                                                         L handle) const {
	for (auto it = pqset.begin(); it != pqset.end(); ++it) {
		if (!condition()) {
			break;
		}
		ElemType elem = *it;
		handle(elem.first, elem.second);
	}
}

} /* namespace Aux */
#endif /* PRIORITYQUEUE_H_ */
