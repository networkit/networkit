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

	// TODO: change this to vector
	T* array; //!< array of size (n+1).  array[0] is not a valid entry, since node indices are 1-based
	T defaultValue; //!< default value
	T nullValue; //!< denotes absence of a value
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
	T& operator[](const I& index);


	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	const T& operator[](const I& index) const;


	/**
	 * Get the number of 1-based entries in this map.
	 */
	inline int64_t numberOfEntries() const;


	/**
	 * Check whether map contains a valid entry for index.
	 */
	inline bool contains(I index) const;


	/**
	 * quick & dirty debug print
	 * TODO: replace with operator<<
	 */
	void print();



};

} /* namespace EnsembleClustering */

template<typename I, typename T>
inline EnsembleClustering::IndexMap<I, T>::IndexMap(int64_t n) {
	this->n = n;
	this->defaultValue = 0;
	this->nullValue = 0;
	this->array = new T[n+1];
	for (int64_t i = 1; i <= n; ++i) {
		this->array[i] = this->nullValue;
	}
}

template<typename I, typename T>
inline EnsembleClustering::IndexMap<I, T>::IndexMap(int64_t n, T defaultValue) {
	this->n = n;
	this->defaultValue = defaultValue;
	this->nullValue = 0;
	this->array = new T[n+1];
	this->array[0] = 0;
	for (int64_t i = 1; i <= n; ++i) {
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

template<typename I, typename T>
inline int64_t EnsembleClustering::IndexMap<I, T>::numberOfEntries() const {
	// assert (this->n == (this->array.size() - 1));
	return this->n;
}

template<typename I, typename T>
inline bool EnsembleClustering::IndexMap<I, T>::contains(I index) const {
	bool cont = (this->array[index] != this->nullValue);
	return cont;
}

template<typename I, typename T>
inline void EnsembleClustering::IndexMap<I, T>::print() {
	std::cout << "{";
	for (int64_t i = 0; i <= this->n; ++i) {
		std::cout << i << ":" << this->array[i] << ", ";
	}
	std::cout << "}" << std::endl;

}

/*** Implementation ***/


#endif /* INDEXMAP_H_ */
