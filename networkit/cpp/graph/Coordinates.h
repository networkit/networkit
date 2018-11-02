/*
 * Coordinates.h
 *
 *  Created on: 27.01.2014
 *      Author: Henning
 */

#ifndef COORDINATES_H_
#define COORDINATES_H_

#include "../viz/Point.h"

#include <vector>


namespace NetworKit {

/**
 * DEPRECATED: A specialized container for node coordinates is no longer consistent with NetworKit's approach
 * to node attributes. Use standard library containers instead.
 *
 * @ingroup graph
 * The Coordinates class represents a vector of Points of elements of type @a T. The class has convenience
 * functions to get the minimum and maximum coordinate.
 */
template<class T>
class Coordinates {
protected:
	std::vector<Point<T> > data;


public:
	/**
	 * Allocates space for @a numNodes entries.
	 */
	void init(count numNodes) {
		data.resize(numNodes);
	}

	/**
	 * Sets entry at index @a v to value @a value.
	 */
	void setCoordinate(index v, const Point<T>& value) {
		data[v] = value;
	}

	/**
	 * @return Entry at index @a v.
	 */
	Point<T>& getCoordinate(index v) {
		return data[v];
	}

	/**
	 * @return Minimum value of all coordinates with respect to dimension @a dim.
	 */
	T minCoordinate(count dim) {
		T value = data[0][dim];
		for (index i = 1; i < data.size(); ++i) {
			T temp = data[i][dim];
			if (temp < value) {
				value = temp;
			}
		}
		return value;
	}

	/**
	 * @return Maximum value of all coordinates with respect to dimension @a dim.
	 */
	T maxCoordinate(count dim) {
		T value = data[0][dim];
		for (index i = 1; i < data.size(); ++i) {
			T temp = data[i][dim];
			if (temp > value) {
				value = temp;
			}
		}
		return value;
	}

	/**
	 * Insert coordinates of a new vertex.
	 */
	void addCoordinates(std::vector<T>& values) {
		Point<T> p(values);
		data.push_back(p);
	}

	virtual ~Coordinates() {}
};


} /* namespace NetworKit */
#endif /* COORDINATES_H_ */
