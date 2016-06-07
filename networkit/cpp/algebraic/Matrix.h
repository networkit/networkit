/*
 * Matrix.h
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "../graph/Graph.h"
#include "Vector.h"
#include "SparseAccumulator.h"
#include "AlgebraicGlobals.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * The matrix class represents a matrix which is optimized for sparse matrices.
 */
class Matrix {
protected:
	Graph graph;

	count nRows;
	count nCols;

	double zero;

public:
	/** Default constructor */
	Matrix();

	/**
	 * Constructs the Matrix with size @a dimension x @a dimension.
	 * @param dimension Defines how many rows and columns this matrix has.
	 * @param zero The zero element (default is 0.0).
	 */
	Matrix(const count dimension, const double zero = 0.0);


	/**
	 * Constructs the Matrix with size @a nRows x @a nCols.
	 * @param nRows Number of rows.
	 * @param nCols Number of columns.
	 * @param zero The zero element (default is 0.0).
	 */
	Matrix(const count nRows, const count nCols, const double zero = 0.0);

	/**
	 * Constructs the @a dimension x @a dimension Matrix from the elements at position @a positions with values @values.
	 * @param dimension Defines how many rows and columns this matrix has.
	 * @param triplets The nonzero elements.
	 * @param zero The zero element (default is 0.0).
	 */
	Matrix(const count dimension, const std::vector<Triplet>& triplets, const double zero = 0.0);

	/**
	 * Constructs the @a nRows x @a nCols Matrix from the elements at position @a positions with values @values.
	 * @param nRows Defines how many rows this matrix has.
	 * @param nCols Defines how many columns this matrix has.
	 * @param triplets The nonzero elements.
	 * @param zero The zero element (default is 0.0).
	 */
	Matrix(const count nRows, const count nCols, const std::vector<Triplet>& triplets, const double zero = 0.0);

	/** Default copy constructor */
	Matrix(const Matrix &other) = default;

	/** Default move constructor */
	Matrix(Matrix &&other) = default;

	/** Default destructor */
	virtual ~Matrix() = default;

	/** Default move assignment operator */
	Matrix& operator=(Matrix &&other) = default;

	/** Default copy assignment operator */
	Matrix& operator=(const Matrix &other) = default;

	bool operator==(const Matrix& other) const {
		bool graphsEqual = graph.numberOfNodes() == other.graph.numberOfNodes() && graph.numberOfEdges() == other.graph.numberOfEdges();
		if (graphsEqual) {
			graph.forEdges([&](node u, node v, edgeweight w) {
				if (w != other.graph.weight(u, v)) {
					graphsEqual = false;
					return;
				}
			});
		}

		return graphsEqual && nRows == other.nRows && nCols == other.nCols && zero == other.zero;
	}

	bool operator!=(const Matrix& other) const {
		return !((*this) == other);
	}

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
	 */
	void setValue(const index i, const index j, const double value);

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
	Matrix operator+(const Matrix &other) const;

	/**
	 * Adds @a other to this matrix.
	 * @return Reference to this matrix.
	 */
	Matrix& operator+=(const Matrix &other);

	/**
	 * Subtracts @a other from this matrix and returns the result.
	 * @return The difference of this matrix and @a other.
	 *
	 */
	Matrix operator-(const Matrix &other) const;

	/**
	 * Subtracts @a other from this matrix.
	 * @return Reference to this matrix.
	 */
	Matrix& operator-=(const Matrix &other);

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar and returns the result.
	 * @return The result of multiplying this matrix with @a scalar.
	 */
	Matrix operator*(const double scalar) const;

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar.
	 * @return Reference to this matrix.
	 */
	Matrix& operator*=(const double scalar);

	/**
	 * Multiplies this matrix with @a vector and returns the result.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;

	/**
	 * Multiplies this matrix with @a other and returns the result in a new matrix.
	 * @return The result of multiplying this matrix with @a other.
	 */
	Matrix operator*(const Matrix &other) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor and returns the result in a new matrix.
	 * @return The result of dividing this matrix by @a divisor.
	 */
	Matrix operator/(const double divisor) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor.
	 * @return Reference to this matrix.
	 */
	Matrix& operator/=(const double divisor);

	/**
	 * Computes A^T * B.
	 * @param A
	 * @param B
	 */
	static Matrix mTmMultiply(const Matrix &A, const Matrix &B);

	/**
	 * Computes A * B^T.
	 * @param A
	 * @param B
	 */
	static Matrix mmTMultiply(const Matrix &A, const Matrix &B);

	/**
	 * Computes matrix^T * vector
	 * @param matrix
	 * @param vector
	 */
	static Vector mTvMultiply(const Matrix &matrix, const Vector &vector);

	/**
	 * Transposes this matrix and returns it.
	 */
	Matrix transpose() const;

	/**
	 * Extracts a matrix with rows and columns specified by @a rowIndices and @a columnIndices from this matrix.
	 * The order of rows and columns is equal to the order in @a rowIndices and @a columnIndices. It is also
	 * possible to specify a row or column more than once to get duplicates.
	 * @param rowIndices
	 * @param columnIndices
	 */
	Matrix extract(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices) const;

	/**
	 * Assign the contents of the matrix @a source to this matrix at rows and columns specified by @a rowIndices and
	 * @a columnIndices. That is, entry (i,j) of @a source is assigned to entry (rowIndices[i], columnIndices[j]) of
	 * this matrix. Note that the dimensions of @rowIndices and @a columnIndices must coincide with the number of rows
	 * and columns of @a source.
	 * @param rowIndices
	 * @param columnIndices
	 * @param source
	 */
	void assign(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices, const Matrix& source);

	/**
	 * Applies the unary function @a unaryElementFunction to each value in the matrix. Note that it must hold that f(0) = 0.
	 * @param unaryElementFunction
	 */
	template<typename F>
	void apply(const F unaryElementFunction);

	/**
	 * Returns the (weighted) adjacency matrix of the (weighted) Graph @a graph.
	 * @param graph
	 */
	static Matrix adjacencyMatrix(const Graph& graph, double zero = 0.0);

	/**
	 * Creates a diagonal matrix with dimension equal to the dimension of the Vector @a diagonalElements. The values on
	 * the diagonal are the ones stored in @a diagonalElements (i.e. D(i,i) = diagonalElements[i]).
	 * @param diagonalElements
	 */
	static Matrix diagonalMatrix(const Vector& diagonalElements, double zero = 0.0);

	/**
	 * Returns the (weighted) incidence matrix of the (weighted) Graph @a graph.
	 * @param graph
	 */
	static Matrix incidenceMatrix(const Graph& graph, double zero = 0.0);

	/**
	 * Returns the (weighted) Laplacian matrix of the (weighteD) Graph @a graph.
	 * @param graph
	 */
	static Matrix laplacianMatrix(const Graph& graph,double zero = 0.0);

	/**
	 * Returns the (weighted) normalized Laplacian matrix of the (weighted) Graph @a graph
	 * @param graph
	 */
	static Matrix normalizedLaplacianMatrix(const Graph& graph, double zero = 0.0);

	/**
	 * Iterate over all non-zero elements of row @a row in the matrix and call handler(index row, index column, double value)
	 */
	template<typename L> void forNonZeroElementsInRow(index row, L handle) const;

	template<typename L> void forElementsInRow(index i, L handle) const;

	/**
	 * Iterate over all non-zero elements of the matrix in row order and call handler (lambda closure).
	 */
	template<typename L> void forNonZeroElementsInRowOrder(L handle) const;

	/**
	 * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
	 */
	template<typename L> void parallelForNonZeroElementsInRowOrder(L handle) const;
};


} /* namespace NetworKit */

template<typename F>
void NetworKit::Matrix::apply(const F unaryElementFunction) {
	forNonZeroElementsInRowOrder([&](index i, index j, double value) {
		setValue(i,j, unaryElementFunction(value));
	});
}

template<typename L>
inline void NetworKit::Matrix::forNonZeroElementsInRow(index row, L handle) const {
	graph.forEdgesOf(row, [&](index j, edgeweight weight){
		handle(j, weight);
	});
}

template<typename L>
inline void NetworKit::Matrix::forElementsInRow(index i, L handle) const {
	Vector rowVector = row(i);
	index j = 0;
	rowVector.forElements([&](double value) {
		handle(j++, value);
	});
}

template<typename L>
inline void NetworKit::Matrix::forNonZeroElementsInRowOrder(L handle) const {
	for (index i = 0; i < nRows; ++i) {
		graph.forEdgesOf(i, [&](index j, edgeweight weight){
			handle(i, j, weight);
		});
	}
}

template<typename L>
inline void NetworKit::Matrix::parallelForNonZeroElementsInRowOrder(L handle) const {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		graph.forEdgesOf(i, [&](index j, edgeweight weight){
			handle(i, j, weight);
		});
	}
}

#endif /* MATRIX_H_ */
