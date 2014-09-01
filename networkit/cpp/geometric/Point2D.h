/*
 * Point2D.h
 *
 *  Created on: 24.07.2014
 *      Author: moritzl
 */

#ifndef POINT2D_H_
#define POINT2D_H_

#include <vector>
#include <cinttypes>
#include <cassert>
#include <cmath>
#include <cstdint>

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity

//template<class T> class Point;
//
//template<class T>
//std::ostream& operator <<(std::ostream& out, Point<T>& point);



/**
 * @ingroup viz
 * Points in any dimension of templated type.
 */
template<class T>
class Point2D {
protected:
	T x;
	T y;

public:
	Point2D()	{
		x = 0;
		y = 0;
	}
	Point2D(T x, T y) {
		this->x = x;
		this->y = y;
	}

	virtual ~Point2D() {}

	count getDimensions() { return 2; }

	T distance(const Point2D<T>& p) const;
	T squaredDistance(const Point2D<T>& p) const;

	Point2D& operator+=(const Point2D<T>& p);
	Point2D& operator-=(const Point2D<T>& p);
	Point2D& scale(const T factor);

	Point2D operator-(const Point2D<T>& other);
	Point2D operator+(const Point2D<T>& other);

	T length() const;
	T squaredLength() const;

	T& operator[](const index i);
};

template<class T>
T Point2D<T>::length() const {
	return sqrt(x*x+y*y);
}

template<class T>
T Point2D<T>::squaredLength() const {
	return x*x + y*y;
}

template<class T>
T Point2D<T>::squaredDistance(const Point2D<T>& p) const {
	T diffx = p.x - x;
	T diffy = p.y - y;
	return diffx*diffx + diffy*diffy;
}

template<class T>
T Point2D<T>::distance(const Point2D<T>& p) const {
	return sqrt(squaredDistance(p));
}

template<class T>
Point2D<T>& Point2D<T>::operator+=(const Point2D<T>& p) {
	this->x += p.x;
	this->y += p.y;
	return *this;
}

template<class T>
Point2D<T>& Point2D<T>::operator-=(const Point2D<T>& p) {
	this->x -= p.x;
	this->y -= p.y;
	return *this;
}

template<class T>
Point2D<T> Point2D<T>::operator-(const Point2D<T>& other) {
	return Point2D(x - other.x, y - other.y);
}

template<class T>
Point2D<T> Point2D<T>::operator+(const Point2D<T>& other) {
	return Point2D(x + other.x, y + other.y);
}


template<class T>
Point2D<T>& Point2D<T>::scale(const T factor) {
	x *= factor;
	y *= factor;
	return *this;
}

template<class T>
inline T& Point2D<T>::operator [](index i) {
	assert(i >= 0 && i <2);
	if (i == 0) return x;
	else return y;
}

} /* namespace NetworKit */
#endif /* POINT2D_H_ */
