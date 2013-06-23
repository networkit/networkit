/*
 * NodeMap.h
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NODEMAP_H_
#define NODEMAP_H_

#include "Graph.h"
#include "../base/IndexMap.h"

namespace NetworKit {

// forward declarations for specialization
template<typename T> class NodeMap;
template<typename T> std::ostream& operator<<(std::ostream& os, const NodeMap<T>& m);


template <typename T> class NodeMap : public IndexMap<node, T> {


public:

	/**
	 * Construct a new IndexMap which holds n entries .
	 *
	 * @param[in]	n				number of entries
	 * @param[in]	defaultValue	all entries are initialized to this value
	 */
	NodeMap(count n, T defaultValue);

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


	// TODO: string representation of NodeMap

	//template<typename T>
	friend std::ostream& operator<< <> (std::ostream& os, const NodeMap<T>& m);
};

} /* namespace NetworKit */


/*** Implementation ***/



template<typename T> inline NetworKit::NodeMap<T>::NodeMap(count n, T defaultValue = 0) :
		IndexMap<node, T>(n, defaultValue) {
	TRACE("NodeMap initialized with n = " << this->n);
}

template<typename T> inline NetworKit::NodeMap<T>::~NodeMap() {
}

template<typename T> inline T& NetworKit::NodeMap<T>::operator [](const node& u) {
	return this->data[u];
}

template<typename T> inline const T& NetworKit::NodeMap<T>::operator [](
		const node& u) const {
	return this->data[u];
}


// FIXME: linker errors when calling operator<<
template<typename T>
std::ostream& operator <<(std::ostream& os, const NetworKit::NodeMap<T>& m) {
	os << "{ ";
	for (NetworKit::node v = 1; v <= m.n; ++v) {
		T val = m.array[v];
		os << v ;
		os << ": ";
		os << val;
		os << ", ";
	}
	return os << "}";
}

#endif /* NODEMAP_H_ */
