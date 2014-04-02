/*
 * Vector.h
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <cmath>

// forward declaration of Matrix class
class Matrix;

class Vector {
private:
	std::vector<double> values;
	bool transposed;

public:
	Vector();
	Vector(const uint64_t dimension, const double initialValue, const bool transpose = false);
	Vector(const std::vector<double> &values, const bool transpose = false);
	Vector(const std::initializer_list<double> &list);
	Vector(const Vector &other, const bool transpose = false);
	virtual ~Vector();

	/**
	 * @return dimension of vector
	 */
	inline uint64_t getDimension() const {
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
	 * @return Reference to the element at index @a idx.
	 */
	inline double& operator()(const unsigned int &idx) {
		if (idx >= values.size()) {
			throw std::out_of_range("index out of range");
		} else {
			return values[idx];
		}
	}

	/**
	 * @return Constant reference to the element at index @a idx.
	 */
	inline const double& operator()(const unsigned int &idx) const {
		if (idx >= values.size()) {
			throw std::runtime_error("index out of range");
		} else {
			return values[idx];
		}
	}

	/**
	 * Compares this vector and @a other element-wise.
	 * @return true, if this vector is element-wise equal to @a other, otherwise false.
	 */
	bool operator==(const Vector &other) const;

	/**
	 * Compares this vector and @a other element-wise.
	 * @return true, if this vector is element-wise unequal to @a other, otherwise false.
	 */
	bool operator!=(const Vector &other) const;

	/**
	 * Computes the outer product of this vector and @a other.
	 * Note that the dimensions must match and that @a other must be transposed.
	 * @return The result of the outer product.
	 */
//	Matrix outerProduct(const Vector &other) const;

	/**
	 * Computes the outer product of vectors @a v1 and @a v2.
	 * Note that the dimensions must match and that @a v2 must be transposed.
	 * @return The result of the outer product.
	 */
	static Matrix outerProduct(const Vector &v1, const Vector &v2);

	/**
	 * Computes the inner product (dot product) of the vectors @a v1 and @a v2.
	 * @return The result of the inner product.
	 */
	static double innerProduct(const Vector &v1, const Vector &v2);

	/**
	 * Computes the inner product (dot product) of this vector and @a other.
	 * Note that the dimensions must match and that this vector must be transposed.
	 * @return The result of the inner product.
	 */
	double operator*(const Vector &other) const;

	/**
	 * Multiplies this vector with a scalar specified in @a scalar and returns the result.
	 * @return The result of multiplying this vector with @a scalar.
	 */
	Vector operator*(const double &scalar) const;

	/**
	 * Multiplies this vector with a scalar specified in @a scalar.
	 * @return Reference to this vector.
	 */
	Vector& operator*=(const double &scalar);

	/**
	 * Adds this vector to @a other and returns the result.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return The sum of this vector and @a other.
	 */
	Vector operator+(const Vector &other) const;

	/**
	 * Adds @a other to this vector.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return Reference to this vector.
	 */
	Vector& operator+=(const Vector &other);

	/**
	 * Subtracts @a other from this vector and returns the result.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return The difference of this vector and @a other.
	 *
	 */
	Vector operator-(const Vector &other) const;

	/**
	 * Subtracts @a other from this vector.
	 * Note that the dimensions of the vectors have to be the same.
	 * @return Reference to this vector.
	 */
	Vector& operator-=(const Vector &other);

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
	template<typename L> void parallelForElements(L handle); //TODO: Is this iterator useful for anything?

	/**
	 * Iterate in parallel over all elements of the vector and call handler (lambda closure).
	 */
	template<typename L> void parallelForElements(L handle) const; //TODO: Is this iterator useful for anything?

};

/**
 * Multiplies the vector @a v with a scalar specified in @a scalar and returns the result.
 * @return The result of multiplying this vector with @a scalar.
 */
inline Vector operator*(const double &scalar, const Vector &v) {
	return v.operator*(scalar);
}

template<typename L>
inline void Vector::forElements(L handle) {
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(values[i]);
	}
}

template<typename L>
inline void Vector::forElements(L handle) const {
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(values[i]);
	}
}

template<typename L>
inline void Vector::parallelForElements(L handle) {
#pragma omp parallel for
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(i, values[i]);
	}
}

template<typename L>
inline void Vector::parallelForElements(L handle) const {
#pragma omp parallel for
	for (uint64_t i = 0; i < getDimension(); i++) {
		handle(i, values[i]);
	}
}

#endif /* VECTOR_H_ */
