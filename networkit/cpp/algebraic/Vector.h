/*
 * Vector.h
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <vector>
#include "../Globals.h"
#include <cassert>

namespace NetworKit {

// forward declaration of Matrix class
class Matrix;

/**
 * @ingroup algebraic
 * The Vector class represents a basic vector with double coefficients.
 */
class Vector {
private:
	std::vector<double> values;
	bool transposed;

public:
	/** Default constructor */
	Vector();

	/**
	 * Constructs the Vector with @a dimension elements with value @a initialValue.
	 * @param dimension The dimension of this vector.
	 * @param initialValue All coefficients will be initialized to @a initialValue.
	 * @param transpose Indicates whether this vector is transposed (row vector) or not (column vector).
	 */
	Vector(const count dimension, const double initialValue = 0, const bool transpose = false);

	/**
	 * Constructs the Vector with the contents of @a values.
	 * @param values The values of this Vector.
	 * @param transpose Indicates whether this vector is transposed (row vector) or not (column vector).
	 */
	Vector(const std::vector<double> &values, const bool transpose = false);

	/**
	 * Constructs the Vector from the contents of the initializer list @a list.
	 * @param list The initializer list.
	 */
	Vector(const std::initializer_list<double> &list);

	/** Default copy constructor */
	Vector(const Vector &other) = default;

	/** Default move constructor */
	Vector(Vector &&other) = default;

	/** Default destructor */
	virtual ~Vector() = default;

	/** Default copy assignment operator */
	Vector& operator=(const Vector &other) = default;

	/** Default move assignment operator */
	Vector& operator=(Vector &&other) = default;


	/**
	 * @return dimension of vector
	 */
	inline count getDimension() const {
		return values.size();
	}

	/**
	 * A transposed vector is a row vector.
	 * @return True, if this vector is transposed, otherwise false.
	 */
	bool isTransposed() const;

	/**
	 * @return Transposed copy of this vector.
	 */
	Vector transpose() const;


	/**
	 * Calculates and returns the Euclidean length of this vector
	 * @return The Euclidean length of this vector.
	 */
	double length() const;


	/**
	 * Calculates and returns the arithmetic mean of this vector
	 * @return The arithmetic mean of this vector.
	 */
	double mean() const;

	/**
	 * Returns a reference to the element at index @a idx without checking the range of this vector.
	 * @param idx The index of the element.
	 * @return Reference to the element at index @a idx.
	 */
	inline double& operator[](const index idx) {
		assert(idx < values.size());
		return values[idx];
	}

	/**
	 * Returns a constant reference to the element at index @a idx without checking the range of this vector.
	 * @a idx The index of the element.
	 * @return Constant reference to the element at index @a idx.
	 */
	inline const double& operator[](const index idx) const {
		assert(idx < values.size());
		return values[idx];
	}

	/**
	 * Returns a reference to the element at index @a idx. If @a idx is not a valid index an exception is thrown.
	 * @param idx The index of the element.
	 * @return Reference to the element at index @a idx.
	 */
	double &at(const index idx) {
		if (idx >= values.size()) {
			throw std::runtime_error("index out of range");
		} else {
			return values[idx];
		}
	}

	/**
	 * Compares this vector and @a other element-wise.
	 * @return True, if this vector is element-wise equal to @a other, otherwise false.
	 */
	bool operator==(const Vector &other) const;

	/**
	 * Compares this vector and @a other element-wise.
	 * @return True, if this vector is element-wise unequal to @a other, otherwise false.
	 */
	bool operator!=(const Vector &other) const;


	/**
	 * Computes the outer product of @a v1 and @a v2.
	 * @param v1 First Vector.
	 * @param v2 Second Vector.
	 * @return The resulting matrix from the outer product.
	 */
	static Matrix outerProduct(const Vector &v1, const Vector &v2);

	/**
	 * Computes the inner product (dot product) of the vectors @a v1 and @a v2.
	 * @return The result of the inner product.
	 */
	static double innerProduct(const Vector &v1, const Vector &v2);

	/**
	 * Computes the inner product (dot product) of this vector and @a other.
	 * @return The result of the inner product.
	 */
	double operator*(const Vector &other) const;

	/**
	 * Multiplies this vector with @a matrix and returns the result.
	 * @return The result of multiplying this vector with @a matrix.
	 */
	Vector operator*(const Matrix &matrix) const;

	/**
	 * Multiplies this vector with a scalar specified in @a scalar and returns the result in a new vector.
	 * @return The result of multiplying this vector with @a scalar.
	 */
	Vector operator*(const double &scalar) const;

	/**
	 * Multiplies this vector with a scalar specified in @a scalar.
	 * @return Reference to this vector.
	 */
	Vector& operator*=(const double &scalar);


	/**
	 * Divides this vector by a divisor specified in @a divisor and returns the result in a new vector.
	 * @return The result of dividing this vector by @a divisor.
	 */
	Vector operator/(const double &divisor) const;

	/**
	 * Divides this vector by a divisor specified in @a divisor.
	 * @return Reference to this vector.
	 */
	Vector& operator/=(const double &divisor);

	/**
	 * Adds this vector to @a other and returns the result.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return The sum of this vector and @a other.
	 */
	Vector operator+(const Vector &other) const;

	/**
	 * Adds @a value to each element of this vector and returns the result.
	 */
	Vector operator+(const double value) const;

	/**
	 * Adds @a other to this vector.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return Reference to this vector.
	 */
	Vector& operator+=(const Vector &other);

	/**
	 * Adds @a value to each element of this vector.
	 */
	Vector& operator+=(const double value);

	/**
	 * Subtracts @a other from this vector and returns the result.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return The difference of this vector and @a other.
	 *
	 */
	Vector operator-(const Vector &other) const;

	/**
	 * Subtracts @a value from each element of this vector and returns the result.
	 */
	Vector operator-(const double value) const;

	/**
	 * Subtracts @a other from this vector.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return Reference to this vector.
	 */
	Vector& operator-=(const Vector &other);

	/**
	 * Subtracts @a value from each element of this vector.
	 */
	Vector& operator-=(const double value);

	/**
	 * Iterate over all elements of the vector and call handler (lambda closure).
	 */
	template<typename L> void forElements(L handle);

	/**
	 * Iterate over all elements of the vector and call handler (lambda closure).
	 */
	template<typename L> void forElements(L handle) const;

	/**
	 * Iterate in parallel over all elements of the vector and call handler (lambda closure).
	 *
	 */
	template<typename L> void parallelForElements(L handle);

	/**
	 * Iterate in parallel over all elements of the vector and call handler (lambda closure).
	 */
	template<typename L> void parallelForElements(L handle) const;

};

} /* namespace NetworKit */

/**
 * Multiplies the vector @a v with a scalar specified in @a scalar and returns the result.
 * @return The result of multiplying this vector with @a scalar.
 */
inline NetworKit::Vector operator*(const double &scalar, const NetworKit::Vector &v) {
	return v.operator*(scalar);
}

template<typename L>
inline void NetworKit::Vector::forElements(L handle) {
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(values[i]);
	}
}

template<typename L>
inline void NetworKit::Vector::forElements(L handle) const {
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(values[i]);
	}
}

template<typename L>
inline void NetworKit::Vector::parallelForElements(L handle) {
#pragma omp parallel for
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(i, values[i]);
	}
}

template<typename L>
inline void NetworKit::Vector::parallelForElements(L handle) const {
#pragma omp parallel for
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(i, values[i]);
	}
}

#endif /* VECTOR_H_ */
