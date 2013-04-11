/*
 * Point.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef POINT_H_
#define POINT_H_

#include <vector>
#include <cinttypes>

namespace NetworKit {

typedef int64_t index; // more expressive name for an index into an array
typedef int64_t count; // more expressive name for an integer quantity


/**
 * Points in any dimension of templated type.
 */
template<class T>
class Point {
protected:
	std::vector<T> data;

public:
	Point(std::vector<T>& values): data(values) {}
	virtual ~Point() {}

	count getDimensions() { return data.size(); }
	T getValue(index dim) { return data.at(dim); }
	void setValue(index dim, T value) { data.at(dim) = value; }
};

} /* namespace NetworKit */
#endif /* POINT_H_ */
