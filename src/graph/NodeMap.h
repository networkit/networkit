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

namespace EnsembleClustering {

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

} /* namespace EnsembleClustering */


/*** Implementation ***/



template<typename T> inline EnsembleClustering::NodeMap<T>::NodeMap(count n, T defaultValue = -1) :
		IndexMap<node, T>(n, defaultValue) {
}

template<typename T> inline EnsembleClustering::NodeMap<T>::~NodeMap() {
}

template<typename T> inline T& EnsembleClustering::NodeMap<T>::operator [](const node& u) {
	return this->data[u];
}

template<typename T> inline const T& EnsembleClustering::NodeMap<T>::operator [](
		const node& u) const {
	return this->data[u];
}


// FIXME: linker errors when calling operator<<
template<typename T>
std::ostream& operator <<(std::ostream& os, const EnsembleClustering::NodeMap<T>& m) {
	os << "{ ";
	for (EnsembleClustering::node v = 1; v <= m.n; ++v) {
		T val = m.array[v];
		os << v ;
		os << ": ";
		os << val;
		os << ", ";
	}
	return os << "}";
}

#endif /* NODEMAP_H_ */
