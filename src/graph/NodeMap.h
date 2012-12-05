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

	T array[];
	T defaultEntry;

public:

	/**
	 * Construct a node map which holds n entries.
	 */
	NodeMap(int n, T defaultEntry = NULL);

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

template<class T> inline EnsembleClustering::NodeMap<T>::NodeMap(int n, T defaultEntry) {
	this->defaultEntry = defaultEntry;
	this->array = {};
	if (defaultEntry != NULL) {
		// initialize all entries to default
		for (int i = 0; i < n; ++i) {
			this->array[i] = defaultEntry;
		}
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
