/*
 * IndexMap.h
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#ifndef INDEXMAP_H_
#define INDEXMAP_H_

namespace EnsembleClustering {

/**
 * An IndexMap implements a 1-based mapping from an integer index type to an arbitray value type.
 *
 */
template <typename I, typename T> class IndexMap {


protected:

	T* array; //!< array of size (n+1).  array[0] is not a valid entry, since node indices are 1-based
	T defaultValue;
	int64_t n;	//<! number of indices

public:

	IndexMap(int64_t n);

	/**
	 * Construct a new IndexMap which holds n entries .
	 *
	 * @param[in]	defaultValue	all entries are initialized to this value
	 */
	IndexMap(int64_t n, T defaultValue);

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

};

} /* namespace EnsembleClustering */

template<typename I, typename T>
inline EnsembleClustering::IndexMap<I, T>::IndexMap(int64_t n) {
	this->n = n;
	this->defaultValue = defaultValue;
	this->array = new T[n+1];
}

template<typename I, typename T>
inline EnsembleClustering::IndexMap<I, T>::IndexMap(int64_t n, T defaultValue) {
	this->n = n;
	this->defaultValue = defaultValue;
	this->array = new T[n+1];
	for (int64_t i = 1; i < n+1; ++i) {
		this->array[i] = defaultValue;
	}
}

template<typename I, typename T>
inline EnsembleClustering::IndexMap<I, T>::~IndexMap() {
	//TODO: destructor stub
}

template<typename I, typename T>
inline T& EnsembleClustering::IndexMap<I, T>::operator [](const I& index) {
	return this->array[index];
}

template<typename I, typename T>
inline const T& EnsembleClustering::IndexMap<I, T>::operator [](const I& index) const {
	return this->array[index];
}

/*** Implementation ***/


#endif /* INDEXMAP_H_ */
