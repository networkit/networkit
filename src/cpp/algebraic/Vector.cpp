/*
 * Vector.cpp
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Vector.h"


Vector::Vector() : values(0), transposed(false) {}

Vector::Vector(const int dimension, const double initialValue, const bool transposed) : values(dimension, initialValue), transposed(transposed) {}

Vector::Vector(const std::vector<double> &values, const bool transposed) : values(values), transposed(transposed) {
}

Vector::~Vector() {
}


double Vector::innerProduct(const Vector &v1, const Vector &v2) {
	if (v1.getDimension() != v2.getDimension()) {
		throw std::runtime_error("dimensions of vectors do not match");
	}

	double result = 0.0;
#pragma omp parallel for reduction(+:result)
	for (int i = 0; i < v1.getDimension(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}

Vector Vector::operator*(const double &scalar) const {
	return Vector(*this) *= scalar;
}

Vector& Vector::operator*=(const double &scalar) {
#pragma omp parallel for
	for (int i = 0; i < getDimension(); i++) {
		values[i] *= scalar;
	}

	return *this;
}

bool Vector::operator==(const Vector &other) const {
	if (getDimension() != other.getDimension()) return false;

	for (int i = 0; i < getDimension(); i++) {
		if (values[i] != other(i)) return false;
	}

	return true;
}

bool Vector::operator!=(const Vector &other) const {
	return !(*this == other);
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
	for (int i = 0; i < getDimension(); i++) {
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
	for (int i = 0; i < getDimension(); i++) {
		values[i] -= other(i);
	}

	return *this;
}



