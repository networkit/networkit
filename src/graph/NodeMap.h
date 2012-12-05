/*
 * NodeMap.h
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#ifndef NODEMAP_H_
#define NODEMAP_H_

#include "Graph.h"

namespace EnsembleClustering {

template <class T> class NodeMap {

protected:

	T* array; //!< array of size (n+1).  array[0] is not a valid entry, since node indices are 1-based
	T defaultValue;

public:

	/**
	 * Construct a node map which holds n entries.
	 */
	NodeMap(int64_t n, T defaultValue);

	virtual ~NodeMap();

	/**
	 *  Index operator.
	 *
	 *  @param[in]	u	a node
	 */
	T& operator[](const node& u);

	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	const T& operator[](const node& u) const;
};

} /* namespace EnsembleClustering */


/*** Implementation ***/

template<class T> inline EnsembleClustering::NodeMap<T>::NodeMap(int64_t n, T defaultValue) {
	this->defaultValue = defaultValue;
	this->array = new T[n+1];
	for (int64_t i = 1; i < n+1; ++i) {
		this->array[i] = defaultValue;
	}
}

template<class T> inline EnsembleClustering::NodeMap<T>::~NodeMap() {
	delete[] this->array;
}

template<class T> inline T& EnsembleClustering::NodeMap<T>::operator [](const node& u) {
	return this->array[u];
}

template<class T> inline const T& EnsembleClustering::NodeMap<T>::operator [](
		const node& u) const {
	return this->array[u];
}

#endif /* NODEMAP_H_ */
