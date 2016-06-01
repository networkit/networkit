/*
 * CSRMatrix.h
 *
 *  Created on: May 6, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef CSRMATRIX_H_
#define CSRMATRIX_H_

#include <vector>
#include "../Globals.h"
#include "AlgebraicGlobals.h"
#include "Vector.h"
#include "../graph/Graph.h"
#include "../algebraic/SparseAccumulator.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * The CSRMatrix class represents a sparse matrix stored in CSR-Format (i.e. compressed sparse row).
 * If speed is important, use this CSRMatrix instead of the Matrix class.
 */
class CSRMatrix {
private:
	std::vector<index> rowIdx;
	std::vector<index> columnIdx;
	std::vector<double> nonZeros;

	count nRows;
	count nCols;
	bool isSorted;
	double zero;

	void quicksort(index left, index right);
	index partition(index left, index right);
	index binarySearchColumns(index left, index right, index j) const;

public:
	/** Default constructor */
	CSRMatrix();

	/**
	 * Constructs the CSRMatrix with size @a dimension x @a dimension.
	 * @param dimension Defines how many rows and columns this matrix has.
	 */
	CSRMatrix(const count dimension, const double zero = 0.0);

	/**
	 * Constructs the CSRMatrix with size @a nRows x @a nCols.
	 * @param nRows Number of rows.
	 * @param nCols Number of columns.
	 */
	CSRMatrix(const count nRows, const count nCols, const double zero = 0.0);

	/**
	 * Constructs the @a dimension x @a dimension Matrix from the elements at position @a positions with values @values.
	 * @param dimension Defines how many rows and columns this matrix has.
	 * @param triplets The nonzero elements.
	 * @param isSorted True, if the triplets are sorted per row. Default is false.
	 */
	CSRMatrix(const count dimension, const std::vector<Triplet>& triplets, const double zero = 0.0, bool isSorted = false);

	/**
	 * Constructs the @a nRows x @a nCols Matrix from the elements at position @a positions with values @values.
	 * @param nRows Defines how many rows this matrix has.
	 * @param nCols Defines how many columns this matrix has.
	 * @param triplets The nonzero elements.
	 * @param isSorted True, if the triplets are sorted per row. Default is false.
	 */
	CSRMatrix(const count nRows, const count nCols, const std::vector<Triplet>& triplets, const double zero = 0.0, bool isSorted = false);

	// TODO: Check if this constructor is used anywhere
	CSRMatrix(const count nRows, const count nCols, const std::vector<std::vector<index>> &columnIdx, const std::vector<std::vector<double>> &values, const double zero = 0.0, bool isSorted = false);

	/**
	 * Constructs the @a nRows x @a nCols Matrix from the elements at position @a positions with values @values.
	 * @param nRows Defines how many rows this matrix has.
	 * @param nCols Defines how many columns this matrix has.
	 * @param rowIdx The rowIdx vector of the CSR format.
	 * @param columnIdx The columnIdx vector of the CSR format.
	 * @param nonZeros The nonZero vector of the CSR format. Should be as long as the @a columnIdx vector.
	 * @param isSorted True, if the triplets are sorted per row. Default is false.
	 */
	CSRMatrix(const count nRows, const count nCols, const std::vector<index>& rowIdx, const std::vector<index>& columnIdx, const std::vector<double>& nonZeros, const double zero = 0.0, bool isSorted = false);

	/** Default copy constructor */
	CSRMatrix (const CSRMatrix &other) = default;

	/** Default move constructor */
	CSRMatrix (CSRMatrix &&other) = default;

	/** Default destructor */
	virtual ~CSRMatrix() = default;

	/** Default move assignment operator */
	CSRMatrix& operator=(CSRMatrix &&other) = default;

	/** Default copy assignment operator */
	CSRMatrix& operator=(const CSRMatrix &other) = default;

	/**
	 * @return Number of rows.
	 */
	inline count numberOfRows() const {
		return nRows;
	}

	/**
	 * @return Number of columns.
	 */
	inline count numberOfColumns() const {
		return nCols;
	}

	/**
	 * Returns the zero element of the matrix.
	 */
	inline double getZero() const {
		return zero;
	}

	/**
	 * @param i The row index.
	 * @return Number of non-zeros in row @a i.
	 */
	count nnzInRow(const index i) const;

	/**
	 * @return Number of non-zeros in this matrix.
	 */
	count nnz() const;

	/**
	 * @return Value at matrix position (i,j).
	 */
	double operator()(const index i, const index j) const;

	/**
	 * Set the matrix at position (@a i, @a j) to @a value.
	 * @note This operation can be linear in the number of non-zeros due to vector element movements
	 */
	void setValue(const index i, const index j, const double value);

	/**
	 * Sorts the column indices in each row for faster access.
	 */
	void sort();

	/**
	 * @return True if the matrix is sorted, otherwise false.
	 */
	bool sorted() const;

	/**
	 * @return Row @a i of this matrix as vector.
	 */
	Vector row(const index i) const;

	/**
	 * @return Column @a j of this matrix as vector.
	 */
	Vector column(const index j) const;

	/**
	 * @return The main diagonal of this matrix.
	 */
	Vector diagonal() const;

	/**
	 * Adds this matrix to @a other and returns the result.
	 * @return The sum of this matrix and @a other.
	 */
	CSRMatrix operator+(const CSRMatrix &other) const;

	/**
	 * Adds @a other to this matrix.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator+=(const CSRMatrix &other);

	/**
	 * Subtracts @a other from this matrix and returns the result.
	 * @return The difference of this matrix and @a other.
	 *
	 */
	CSRMatrix operator-(const CSRMatrix &other) const;

	/**
	 * Subtracts @a other from this matrix.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator-=(const CSRMatrix &other);

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar and returns the result.
	 * @return The result of multiplying this matrix with @a scalar.
	 */
	CSRMatrix operator*(const double &scalar) const;

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator*=(const double &scalar);

	/**
	 * Multiplies this matrix with @a vector and returns the result.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;

	/**
	 * Multiplies this matrix with @a other and returns the result in a new matrix.
	 * @return The result of multiplying this matrix with @a other.
	 */
	CSRMatrix operator*(const CSRMatrix &other) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor and returns the result in a new matrix.
	 * @return The result of dividing this matrix by @a divisor.
	 */
	CSRMatrix operator/(const double &divisor) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator/=(const double &divisor);

	/**
	 * Creates a submatrix of this matrix consisting of the rows specified in @a rows and columns specified in @a columns.
	 * @param rows The row indices referencing the rows to include in the submatrix.
	 * @param columns The column indices referencing the columns to include in the submatrix.
	 * @return The submatrix of this matrix consisting of @a rows and @a columns.
	 */
	CSRMatrix subMatrix(const std::vector<index> &rows, const std::vector<index> &columns) const;

	/**
	 * Computes @a A @a binaryOp @a B on the elements of matrix @a A and matrix @a B.
	 * @param A Sorted CSRMatrix.
	 * @param B Sorted CSRMatrix.
	 * @param binaryOp Function handling (double, double) -> double
	 * @return @a A @a binaryOp @a B.
	 * @note @a A and @a B must have the same dimensions and must be sorted.
	 */
	template<typename L> static CSRMatrix binaryOperator(const CSRMatrix &A, const CSRMatrix &B, L binaryOp);

	/**
	 * Computes @a A^T * @a B.
	 * @param A
	 * @param B
	 * @return @a A^T * @a B.
	 * @note The number of rows of @a A must be equal to the number of rows of @a B.
	 */
	static CSRMatrix mTmMultiply(const CSRMatrix &A, const CSRMatrix &B);

	/**
	 * Computes @a A * @a B^T.
	 * @param A
	 * @param B
	 * @return @a A * @a B^T.
	 * @note The number of columns of @a A must be equal to the number of columns of @a B.
	 */
	static CSRMatrix mmTMultiply(const CSRMatrix &A, const CSRMatrix &B);

	/**
	 * Computes @a matrix^T * @a vector.
	 * @param matrix
	 * @param vector
	 * @return @a matrix^T * @a vector.
	 * @note The number of rows of @a matrix must be equal to the dimension of @a vector.
	 */
	static Vector mTvMultiply(const CSRMatrix &matrix, const Vector &vector);

	/**
	 * Transposes this matrix and returns it.
	 * @return The transposed matrix of this matrix.
	 */
	CSRMatrix transpose() const;

	/**
	 * Compute the (weighted) adjacency matrix of the (weighted) Graph @a graph.
	 * @param graph
	 */
	static CSRMatrix adjacencyMatrix(const Graph& graph);

	/**
	 * Creates a diagonal matrix with dimension equal to the dimension of the Vector @a diagonalElements. The values on
	 * the diagonal are the ones stored in @a diagonalElements (i.e. D(i,i) = diagonalElements[i]).
	 * @param diagonalElements
	 */
	static CSRMatrix diagonalMatrix(const Vector& diagonalElements);

	/**
	 * Returns the (weighted) incidence matrix of the (weighted) Graph @a graph.
	 * @param graph
	 */
	static CSRMatrix incidenceMatrix(const Graph& graph);

	/**
	 * Compute the (weighted) Laplacian of the (weighted) Graph @a graph.
	 * @param graph
	 */
	static CSRMatrix laplacianMatrix(const Graph& graph);

	/**
	 * Returns the (weighted) normalized Laplacian matrix of the (weighted) Graph @a graph
	 * @param graph
	 */
	static CSRMatrix normalizedLaplacianMatrix(const Graph& graph);



	/**
	 * Iterate over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
	 */
	template<typename L> void forNonZeroElementsInRow(index row, L handle) const;

	/**
	 * Iterate in parallel over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
	 */
	template<typename L> void parallelForNonZeroElementsInRow(index row, L handle) const;

	/**
	 * Iterate over all non-zero elements of the matrix in row order and call handler (lambda closure).
	 */
	template<typename L> void forNonZeroElementsInRowOrder(L handle) const;

	/**
	 * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
	 */
	template<typename L> void parallelForNonZeroElementsInRowOrder(L handle) const;
};

template<typename L> inline CSRMatrix NetworKit::CSRMatrix::binaryOperator(const CSRMatrix &A, const CSRMatrix &B, L binaryOp) {
	assert(A.nRows == B.nRows && A.nCols == B.nCols);

	if (A.sorted() && B.sorted()) {
		std::vector<index> rowIdx(A.nRows+1);
		std::vector<std::vector<index>> columns(A.nRows);

		rowIdx[0] = 0;
#pragma omp parallel for
		for (index i = 0; i < A.nRows; ++i) {
			index k = A.rowIdx[i];
			index l = B.rowIdx[i];
			while (k < A.rowIdx[i+1] && l < B.rowIdx[i+1]) {
				if (A.columnIdx[k] < B.columnIdx[l]) {
					columns[i].push_back(A.columnIdx[k]);
					++k;
				} else if (A.columnIdx[k] > B.columnIdx[l]) {
					columns[i].push_back(B.columnIdx[l]);
					++l;
				} else { // A.columnIdx[k] == B.columnIdx[l]
					columns[i].push_back(A.columnIdx[k]);
					++k;
					++l;
				}
				++rowIdx[i+1];
			}

			while (k < A.rowIdx[i+1]) {
				columns[i].push_back(A.columnIdx[k]);
				++k;
				++rowIdx[i+1];
			}

			while (l < B.rowIdx[i+1]) {
				columns[i].push_back(B.columnIdx[l]);
				++l;
				++rowIdx[i+1];
			}
		}


		for (index i = 0; i < A.nRows; ++i) {
			rowIdx[i+1] += rowIdx[i];
		}

		count nnz = rowIdx[A.nRows];
		std::vector<index> columnIdx(nnz);
		std::vector<double> nonZeros(nnz, A.zero);

#pragma omp parallel for
		for (index i = 0; i < A.nRows; ++i) {
			for (index cIdx = rowIdx[i], j = 0; cIdx < rowIdx[i+1]; ++cIdx, ++j) {
				columnIdx[cIdx] = columns[i][j];
			}
			columns[i].clear();
			columns[i].resize(0);
			columns[i].shrink_to_fit();
		}

#pragma omp parallel for
		for (index i = 0; i < A.nRows; ++i) {
			index k = A.rowIdx[i];
			index l = B.rowIdx[i];
			for (index cIdx = rowIdx[i]; cIdx < rowIdx[i+1]; ++cIdx) {
				if (k < A.rowIdx[i+1] && columnIdx[cIdx] == A.columnIdx[k]) {
					nonZeros[cIdx] = A.nonZeros[k];
					++k;
				}

				if (l < B.rowIdx[i+1] && columnIdx[cIdx] == B.columnIdx[l]) {
					nonZeros[cIdx] = binaryOp(nonZeros[cIdx], B.nonZeros[l]);
					++l;
				}
			}
		}

		return CSRMatrix(A.nRows, A.nCols, rowIdx, columnIdx, nonZeros, A.zero, true);
	} else { // A or B not sorted
		std::vector<int64_t> columnPointer(A.nCols, -1);
		std::vector<double> Arow(A.nCols, A.zero);
		std::vector<double> Brow(A.nCols, B.zero);

		std::vector<Triplet> triplets;

		for (index i = 0; i < A.nRows; ++i) {
			index listHead = 0;
			count nnz = 0;

			// search for nonZeros in our own matrix
			for (index k = A.rowIdx[i]; k < A.rowIdx[i+1]; ++k) {
				index j = A.columnIdx[k];
				Arow[j] = A.nonZeros[k];

				columnPointer[j] = listHead;
				listHead = j;
				nnz++;
			}

			// search for nonZeros in the other matrix
			for (index k = B.rowIdx[i]; k < B.rowIdx[i+1]; ++k) {
				index j = B.columnIdx[k];
				Brow[j] = B.nonZeros[k];

				if (columnPointer[j] == -1) { // our own matrix does not have a nonZero entry in column j
					columnPointer[j] = listHead;
					listHead = j;
					nnz++;
				}
			}

			// apply operator on the found nonZeros in A and B
			for (count k = 0; k < nnz; ++k) {
				double value = binaryOp(Arow[listHead], Brow[listHead]);
				if (value != A.zero) {
					triplets.push_back({i,listHead,value});
				}

				index temp = listHead;
				listHead = columnPointer[listHead];

				// reset for next row
				columnPointer[temp] = -1;
				Arow[temp] = A.zero;
				Brow[temp] = B.zero;
			}

			nnz = 0;
		}

		return CSRMatrix(A.numberOfRows(), A.numberOfColumns(), triplets);
	}


}

} /* namespace NetworKit */

template<typename L>
inline void NetworKit::CSRMatrix::forNonZeroElementsInRow(index i, L handle) const {
	for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
		handle(columnIdx[k], nonZeros[k]);
	}
}

template<typename L>
inline void NetworKit::CSRMatrix::parallelForNonZeroElementsInRow(index i, L handle) const {
#pragma omp parallel for
	for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
		handle(columnIdx[k], nonZeros[k]);
	}
}

template<typename L>
inline void NetworKit::CSRMatrix::forNonZeroElementsInRowOrder(L handle) const {
	for (index i = 0; i < nRows; ++i) {
		for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
			handle(i, columnIdx[k], nonZeros[k]);
		}
	}
}

template<typename L>
inline void NetworKit::CSRMatrix::parallelForNonZeroElementsInRowOrder(L handle) const {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
			handle(i, columnIdx[k], nonZeros[k]);
		}
	}
}

#endif /* TESTMATRIX_H_ */
