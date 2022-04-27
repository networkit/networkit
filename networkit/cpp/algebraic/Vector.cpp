/*
 * Vector.cpp
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <networkit/algebraic/Vector.hpp>

#include <networkit/algebraic/DynamicMatrix.hpp>

namespace NetworKit {

Vector::Vector() : values(0), transposed(false) {}

Vector::Vector(const count dimension, const double initialValue, const bool transpose)
    : values(dimension, initialValue), transposed(transpose) {}

Vector::Vector(const std::vector<double> &values, const bool transpose)
    : values(values), transposed(transpose) {}

Vector::Vector(const std::initializer_list<double> &list) : values(list), transposed(false) {}

bool Vector::isTransposed() const {
    return transposed;
}

Vector Vector::transpose() const {
    Vector v(*this);
    v.transposed = !this->transposed;
    return v;
}

double Vector::length() const {
    return std::sqrt(innerProduct(*this, *this));
}

double Vector::mean() const {
    double sum = 0.;
    const auto n = getDimension();
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(+ : sum)
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i)
        sum += values[i];
#else
    sum = std::accumulate(values.begin(), values.end(), 0.);
#endif

    return sum / (double)n;
}

bool Vector::operator==(const Vector &other) const {
    return isTransposed() == other.isTransposed() && values == other.values;
}

bool Vector::operator!=(const Vector &other) const {
    return !(*this == other);
}

double Vector::innerProduct(const Vector &v1, const Vector &v2) {
    assert(v1.getDimension() == v2.getDimension());
    double inner_prod = 0.;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(+ : inner_prod)
    for (omp_index i = 0; i < static_cast<omp_index>(v1.getDimension()); ++i)
        inner_prod += v1[i] * v2[i];
#else
    inner_prod = std::inner_product(v1.values.begin(), v1.values.end(), v2.values.begin(), 0.0);
#endif
    return inner_prod;
}

double Vector::operator*(const Vector &other) const {
    assert(isTransposed()
           && !other.isTransposed()); // vectors must be transposed correctly for inner product
    assert(getDimension() == other.getDimension()); // dimensions of vectors must match

    return innerProduct(*this, other);
}

Vector Vector::operator*(double scalar) const {
    return Vector(*this) *= scalar;
}

Vector &Vector::operator*=(double scalar) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); i++) {
        values[i] *= scalar;
    }

    return *this;
}

Vector Vector::operator/(double divisor) const {
    return Vector(*this) /= divisor;
}

Vector &Vector::operator/=(double divisor) {
    return *this *= 1 / divisor;
}

Vector Vector::operator+(const Vector &other) const {
    return Vector(*this) += other;
}

Vector Vector::operator+(double value) const {
    return Vector(*this) += value;
}

Vector &Vector::operator+=(const Vector &other) {
    assert(isTransposed() == other.isTransposed()); // vectors must be transposed correctly
    assert(getDimension() == other.getDimension()); // dimensions of vectors must match

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); i++) {
        values[i] += other[i];
    }

    return *this;
}

Vector &Vector::operator+=(const double value) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); ++i) {
        values[i] += value;
    }

    return *this;
}

Vector Vector::operator-(const Vector &other) const {
    return Vector(*this) -= other;
}

Vector Vector::operator-(const double value) const {
    return Vector(*this) -= value;
}

Vector &Vector::operator-=(const Vector &other) {
    assert(isTransposed() == other.isTransposed()); // vectors must be transposed correctly
    assert(getDimension() == other.getDimension()); // dimensions of vectors must match

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); i++) {
        values[i] -= other[i];
    }

    return *this;
}

Vector &Vector::operator-=(const double value) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); ++i) {
        values[i] -= value;
    }

    return *this;
}

} /* namespace NetworKit */
