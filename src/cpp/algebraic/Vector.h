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

class Vector {
private:
	std::vector<double> values;

public:
	Vector();
	Vector(const int dimension, const double initialValue);
	Vector(const std::vector<double> &values);
	virtual ~Vector();

	/**
	 * @return dimension of vector
	 */
	inline int getDimension() const {
		return values.size();
	}

	// TODO: inner and outer product

	/**
	 * @return Reference to the element at index @a idx
	 */
	inline double& operator()(const unsigned int &idx) {
		if (idx >= values.size()) {
			throw std::runtime_error("index out of range");
		} else {
			return values[idx];
		}
	}

	/**
	 * @return Reference to the element at index @a idx
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
	for (int i = 0; i < getDimension(); i++) {
		handle(values[i]);
	}
}

template<typename L>
inline void Vector::forElements(L handle) const {
	for (int i = 0; i < getDimension(); i++) {
		handle(values[i]);
	}
}

template<typename L>
inline void Vector::parallelForElements(L handle) {
#pragma omp parallel for
	for (int i = 0; i < getDimension(); i++) {
		handle(i, values[i]);
	}
}

template<typename L>
inline void Vector::parallelForElements(L handle) const {
#pragma omp parallel for
	for (int i = 0; i < getDimension(); i++) {
		handle(i, values[i]);
	}
}

#endif /* VECTOR_H_ */
