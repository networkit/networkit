/*
 * Vector.cpp
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Vector.h"


Vector::Vector() : values(0), transposed(false) {}

Vector::Vector(const uint64_t dimension, const double initialValue, const bool transpose) : values(dimension, initialValue), transposed(transpose) {}

Vector::Vector(const std::vector<double> &values, const bool transpose) : values(values), transposed(transpose) {
}

Vector::Vector(const std::initializer_list<double> &list) : values(list), transposed(false) {
}

Vector::~Vector() {
}

bool Vector::isTransposed() const {
	return transposed;
}

void Vector::transpose() {
	transposed = !transposed;
}

bool Vector::operator==(const Vector &other) const {
	if (getDimension() != other.getDimension()) return false;

	for (uint64_t i = 0; i < getDimension(); i++) {
		if (values[i] != other(i)) return false;
	}

	return true;
}

bool Vector::operator!=(const Vector &other) const {
	return !(*this == other);
}

//Matrix Vector::outerProduct(const Vector &other) const {
//	return outerProduct(*this, other);
//}
//
//Matrix Vector::outerProduct(const Vector &v1, const Vector &v2) {
//	// TODO: Problem: Matrix class only supports symmetric matrices.
//}

double Vector::operator*(const Vector &other) const {
	if (!transposed || other.isTransposed()) {
		throw std::runtime_error("vectors are not transposed correctly for inner product");
	} else if (getDimension() != other.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

	double result = 0.0;
#pragma omp parallel for reduction(+:result)
	for (uint64_t i = 0; i < getDimension(); ++i) {
		result += (*this)(i) * other(i);
	}

	return result;
}

Vector Vector::operator*(const double &scalar) const {
	return Vector(*this) *= scalar;
}

Vector& Vector::operator*=(const double &scalar) {
#pragma omp parallel for
	for (uint64_t i = 0; i < getDimension(); i++) {
		values[i] *= scalar;
	}

	return *this;
}

Vector Vector::operator+(const Vector &other) const {
	if (this->getDimension() != other.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

	return Vector(*this) += other;
}

Vector& Vector::operator+=(const Vector &other) {
	if (getDimension() != other.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

#pragma omp parallel for
	for (uint64_t i = 0; i < getDimension(); i++) {
		values[i] += other(i);
	}

	return *this;
}

Vector Vector::operator-(const Vector &other) const {
	if (getDimension() != other.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

	return Vector(*this) -= other;
}

Vector& Vector::operator-=(const Vector &other) {
	if (getDimension() != other.getDimension()) {
			throw std::runtime_error("dimensions of vectors do not match");
	}

#pragma omp parallel for
	for (uint64_t i = 0; i < getDimension(); i++) {
		values[i] -= other(i);
	}

	return *this;
}



