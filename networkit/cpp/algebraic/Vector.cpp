/*
 * Vector.cpp
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Vector.h"
#include "Matrix.h"

namespace NetworKit {

Vector::Vector() : values(0), transposed(false) {}

Vector::Vector(const uint64_t dimension, const double initialValue, const bool transpose) : values(dimension, initialValue), transposed(transpose) {}

Vector::Vector(const std::vector<double> &values, const bool transpose) : values(values), transposed(transpose) {
}

Vector::Vector(const std::initializer_list<double> &list) : values(list), transposed(false) {
}

Vector::Vector(const Vector &other) : values(other.values), transposed(other.transposed) {
}

bool Vector::isTransposed() const {
	return transposed;
}

Vector Vector::transpose() const {
	Vector v(*this);
	v.transposed = !this->transposed;
	return v;
}

double Vector::length() const {
	return std::sqrt(this->transpose() * (*this));
}

bool Vector::operator==(const Vector &other) const {
	if (getDimension() != other.getDimension() || isTransposed() != other.isTransposed()) return false;

	for (count i = 0; i < getDimension(); i++) {
		if (values[i] != other[i]) return false;
	}

	return true;
}

bool Vector::operator!=(const Vector &other) const {
	return !(*this == other);
}

double Vector::operator*(const Vector &other) const {
	if (!isTransposed() || other.isTransposed()) {
		throw std::runtime_error("vectors are not transposed correctly for inner product");
	} else if (getDimension() != other.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

	double result = 0.0;
#pragma omp parallel for reduction(+:result)
	for (count i = 0; i < getDimension(); ++i) {
		result += values[i] * other[i];
	}

	return result;
}

Vector Vector::operator*(const Matrix &matrix) const {
	if (!isTransposed()) {
		throw std::runtime_error("vector must be of the form 1xn");
	} else if (getDimension() != matrix.numberOfRows()) {
		throw std::runtime_error("dimensions of vector and matrix do not match");
	}

	Vector result(matrix.numberOfColumns(), 0.0, true);
#pragma omp parallel for
	for (count k = 0; k < matrix.numberOfColumns(); ++k) {
		Vector column = matrix.column(k);
		result[k] = (*this) * column;
	}

	return result;
}

Vector Vector::operator*(const double &scalar) const {
	return Vector(*this) *= scalar;
}

Vector& Vector::operator*=(const double &scalar) {
#pragma omp parallel for
	for (count i = 0; i < getDimension(); i++) {
		values[i] *= scalar;
	}

	return *this;
}

Vector Vector::operator/(const double &divisor) const {
	return Vector(*this) /= divisor;
}

Vector& Vector::operator/=(const double &divisor) {
	return *this *= 1 / divisor;
}

Vector Vector::operator+(const Vector &other) const {
	return Vector(*this) += other;
}

Vector& Vector::operator+=(const Vector &other) {
	if (isTransposed() != other.isTransposed()) {
		throw std::runtime_error("vectors are not transformed correctly");
	} else if (getDimension() != other.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

#pragma omp parallel for
	for (count i = 0; i < getDimension(); i++) {
		values[i] += other[i];
	}

	return *this;
}

Vector Vector::operator-(const Vector &other) const {
	return Vector(*this) -= other;
}

Vector& Vector::operator-=(const Vector &other) {
	if (isTransposed() != other.isTransposed()) {
		throw std::runtime_error("vectors are not transformed correctly");
	} else if (getDimension() != other.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

#pragma omp parallel for
	for (count i = 0; i < getDimension(); i++) {
		values[i] -= other[i];
	}

	return *this;
}


} /* namespace NetworKit */

