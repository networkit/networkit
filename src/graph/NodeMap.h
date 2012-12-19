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
	int64_t n;	//<! number of indices

public:

	/**
	 * Construct a node map which holds n entries. Node ids are 1-based.
	 */
	NodeMap(int64_t n);

	/**
	 * Construct a node map which holds n entries .
	 *
	 * @param[in]	defaultValue	all entries are initialized to this value
	 */
	NodeMap(int64_t n, T defaultValue);

	virtual ~NodeMap();

	/**
	 *  Index operator.
	 *
	 *  @param[in]	u	a node
	 */
	inline T& operator[](const node& u);

	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	inline const T& operator[](const node& u) const;


//	FIXME: friend std::ostream& operator <<(std::ostream& os, const NodeMap<T> m);
};

} /* namespace EnsembleClustering */


/*** Implementation ***/


template<class T> inline EnsembleClustering::NodeMap<T>::NodeMap(int64_t n) {
	this->n = n;
	this->array = new T[n+1];
}

template<class T> inline EnsembleClustering::NodeMap<T>::NodeMap(int64_t n, T defaultValue) {
	this->n = n;
	this->defaultValue = defaultValue;
	this->array = new T[n+1];
	for (int64_t i = 1; i <= n; ++i) {
		this->array[i] = defaultValue;
	}
}

template<class T> inline EnsembleClustering::NodeMap<T>::~NodeMap() {
	delete[] array;
}

template<class T> inline T& EnsembleClustering::NodeMap<T>::operator [](const node& u) {
	return this->array[u];
}

template<class T> inline const T& EnsembleClustering::NodeMap<T>::operator [](
		const node& u) const {
	return this->array[u];
}

// FIXME: operator<<
//template<class T>
//inline std::ostream& operator <<(std::ostream& os, const EnsembleClustering::NodeMap<T> m) {
//	os << "{ ";
//	for (EnsembleClustering::node v = 1; v <= m.n; ++v) {
//		T val = m.array[v];
//		os << v ;
//		os << ": ";
//		os << val;
//		os << ", ";
//	}
//	return os << "}";
//}

#endif /* NODEMAP_H_ */
