/*
 * DenseMatrix.cpp
 *
 *  Created on: Nov 25, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "DenseMatrix.h"

namespace NetworKit {

DenseMatrix::DenseMatrix() : nRows(0), nCols(0), entries(std::vector<double>(0)) {
}

DenseMatrix::DenseMatrix(const count nRows, const count nCols, const std::vector<double> &entries) : nRows(nRows), nCols(nCols), entries(entries) {
	assert(entries.size() == nRows * nCols);
}

double DenseMatrix::operator()(const index i, const index j) const {
	return entries[i * numberOfColumns() + j];
}

void DenseMatrix::setValue(const index i, const index j, const double value) {
	entries[i * numberOfColumns() + j] = value;
}

Vector DenseMatrix::row(const index i) const {
	Vector row(numberOfColumns(), 0.0, true);
	index offset = i * numberOfColumns();
#pragma omp parallel for
	for (index j = 0; j < numberOfColumns(); ++j) {
		row[j] = entries[offset + j];
	}

	return row;
}

Vector DenseMatrix::column(const index j) const {
	Vector column(numberOfRows(), 0.0);
#pragma omp parallel for
	for (index i = 0; i < numberOfRows(); ++i) {
		column[i] = entries[i * numberOfColumns() + j];
	}

	return column;
}

Vector DenseMatrix::diagonal() const {
	Vector diagonal(std::min(numberOfRows(), numberOfColumns()), 0.0);
#pragma omp parallel for
	for (index i = 0; i < diagonal.getDimension(); ++i) {
		diagonal[i] = (*this)(i,i);
	}

	return diagonal;
}

DenseMatrix DenseMatrix::operator+(const DenseMatrix &other) const {
	assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
	return DenseMatrix::binaryOperator(*this, other, [](double val1, double val2){return val1 + val2;});
}

DenseMatrix& DenseMatrix::operator+=(const DenseMatrix &other) {
	assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
	*this = DenseMatrix::binaryOperator(*this, other, [](double val1, double val2){return val1 + val2;});
	return *this;
}

DenseMatrix DenseMatrix::operator-(const DenseMatrix &other) const {
	assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
	return DenseMatrix::binaryOperator(*this, other, [](double val1, double val2){return val1 - val2;});
}

DenseMatrix& DenseMatrix::operator-=(const DenseMatrix &other) {
	assert(numberOfRows() == other.numberOfRows() && numberOfColumns() == other.numberOfColumns());
	*this = DenseMatrix::binaryOperator(*this, other, [](double val1, double val2){return val1 + val2;});
	return *this;
}

DenseMatrix DenseMatrix::operator*(const double &scalar) const {
	return DenseMatrix(*this) *= scalar;
}

DenseMatrix& DenseMatrix::operator*=(const double &scalar) {
#pragma omp parallel for
	for (index k = 0; k < entries.size(); ++k) {
		entries[k] *= scalar;
	}

	return *this;
}

Vector DenseMatrix::operator*(const Vector &vector) const {
	assert(!vector.isTransposed());
	assert(numberOfColumns() == vector.getDimension());

	Vector result(numberOfRows(), 0.0);
#pragma omp parallel for
	for (index i = 0; i < numberOfRows(); ++i) {
		index offset = i * numberOfColumns();
		for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
			result[i] += entries[k] * vector[j];
		}
	}

	return result;
}

DenseMatrix DenseMatrix::operator*(const DenseMatrix &other) const {
	assert(numberOfColumns() == other.numberOfRows());
	std::vector<double> resultEntries(numberOfRows() * other.numberOfColumns());

#pragma omp parallel for
	for (index i = 0; i < numberOfRows(); ++i) {
		index offset = i * numberOfRows();
		for (index k = 0; k < numberOfColumns(); ++k) {
			double val_i_k = (*this)(i,k);
			for (index j = 0; j < other.numberOfColumns(); ++j) {
				resultEntries[offset + j] += val_i_k * other(k,j);
			}
		}
	}

	return DenseMatrix(numberOfRows(), other.numberOfColumns(), resultEntries);
}

DenseMatrix DenseMatrix::operator/(const double &divisor) const {
	return DenseMatrix(*this) /= divisor;
}

DenseMatrix& DenseMatrix::operator/=(const double &divisor) {
	return *this *= 1.0 / divisor;
}

void DenseMatrix::LUDecomposition(DenseMatrix &matrix) {
	assert(matrix.numberOfRows() == matrix.numberOfColumns());
	for (index k = 0; k < matrix.numberOfRows()-1; ++k) {
		for (index i = k+1; i < matrix.numberOfRows(); ++i) {
			matrix.setValue(i, k, matrix(i,k) / matrix(k,k));
			for (index j = k+1; j < matrix.numberOfRows(); ++j) {
				matrix.setValue(i, j, matrix(i,j) - (matrix(i,k) * matrix(k,j)));
			}
		}
	}
}

Vector DenseMatrix::LUSolve(const DenseMatrix &LU, const Vector &b) {
	Vector x = b;

	for (index i = 0; i < LU.numberOfRows()-1; ++i) { // forward substitution
		for (index j = i+1; j < LU.numberOfRows(); ++j) {
			x[j] -= x[i] * LU(j,i);
		}
	}

	for (index i = LU.numberOfRows(); i-- > 0;) { // backward substitution
		x[i] /= LU(i,i);
		for (index j = 0; j < i; ++j) {
			x[j] -= x[i] * LU(j,i);
		}
	}

	return x;
}



} /* namespace NetworKit */
