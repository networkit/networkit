/*
 * IndexMap.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef INDEXMAP_H_
#define INDEXMAP_H_

#include <vector>
#include "../auxiliary/Log.h"

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity


/**
 * @DEPRECATED: This class is deprecated. Do not use it in new code.
 *
 * An IndexMap implements a 0-based mapping from an integer index type to an arbitray value type.
 *
 */
template <typename I, typename T> class IndexMap {


protected:

	count n;	//<! number of indices
	std::vector<T> data; //!< array of size (n).
	T defaultValue; //!< default value

public:

//	IndexMap(int64_t n);

	/**
	 * Construct a new IndexMap which holds n entries .
	 *
	 * @param[in]	n				number of entries
	 * @param[in]	defaultValue	all entries are initialized to this value
	 */
	IndexMap(count n, T defaultValue = 0);

	virtual ~IndexMap();


	/**
	 *  Index operator.
	 *
	 *  @param[in]	u	a node
	 */
	inline T& operator[](const I& index);


	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	inline const T& operator[](const I& index) const;


	virtual T at(I index);


	/**
	 * Get the number of 1-based entries in this map.
	 */
	inline count numberOfEntries() const;


	/**
	 * Check whether map contains an entry other than the default.
	 */
	inline bool hasBeenSet(I index) const;


	/**
	 * Number of nodes for which this clsutering can hold an entry.
	 */
	virtual count numberOfNodes() const; // TODO: appropriate here?


	/**
	 * Set all values to one value
	 */
	virtual void setAll(T value);




	/**
	 * quick & dirty debug print
	 * TODO: replace with operator<<
	 */
	void print();



	std::vector<T> getVector();

	/**
	 * Return string representation.
	 */
	// std::string toString();



};



template<typename I, typename T>
inline IndexMap<I, T>::IndexMap(count n, T defaultValue) :
		data(n, defaultValue), defaultValue(defaultValue) {
	TRACE("IndexMap initialized with n = ",n);
	this->n = n;
}

template<typename I, typename T>
inline IndexMap<I, T>::~IndexMap() {

}

template<typename I, typename T>
inline T& IndexMap<I, T>::operator [](const I& index) {
	return this->data[index];
}

template<typename I, typename T>
inline const T& IndexMap<I, T>::operator [](const I& index) const {
	return this->data[index];
}

template<typename I, typename T>
inline count IndexMap<I, T>::numberOfEntries() const {
	// assert (this->n == (this->array.size() - 1));
	return this->n;
}


template<typename I, typename T>
inline bool IndexMap<I, T>::hasBeenSet(I index) const {
	bool cont = (this->data[index] != this->defaultValue);
	return cont;
}

template<typename I, typename T>
inline void IndexMap<I, T>::print() {
	std::cout << "{";
	for (index i = 0; i < this->n; ++i) {
		std::cout << i << ":" << this->data[i] << ", ";
	}
	std::cout << "}" << std::endl;

}

template<typename I, typename T>
inline count IndexMap<I, T>::numberOfNodes() const {
	return this->data.size();	// indices are 0-based
}

template<typename I, typename T>
inline void IndexMap<I, T>::setAll(T value) {
	#pragma omp parallel for if (this->n >= 100000)  // TODO: correct parallelization condition?
	for (index i = 0; i < this->n; ++i) {
		this->data[i] = value;
	}
}

template<typename I, typename T>
inline T IndexMap<I, T>::at(I index) {
	return this->data.at(index);
}


template<typename I, typename T>
inline std::vector<T> NetworKit::IndexMap<I, T>::getVector() {
	return this->data;
}


} /* namespace NetworKit */



#endif /* INDEXMAP_H_ */
