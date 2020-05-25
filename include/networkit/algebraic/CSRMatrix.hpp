/*
 * CSRMatrix.hpp
 *
 *  Created on: May 6, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_CSR_MATRIX_HPP_
#define NETWORKIT_ALGEBRAIC_CSR_MATRIX_HPP_

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/algebraic/SparseAccumulator.hpp>
#include <networkit/auxiliary/Timer.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * The CSRMatrix class represents a sparse matrix stored in CSR-Format (i.e. compressed sparse row).
 * If speed is important, use this CSRMatrix instead of the Matrix class.
 */
class CSRMatrix final {
private:
    std::vector<index> rowIdx;
    std::vector<index> columnIdx;
    std::vector<double> nonZeros;

    count nRows;
    count nCols;
    bool isSorted;
    double zero;

    /**
     * Quicksort algorithm on columnIdx between [@a left, @a right].
     * @param left
     * @param right
     */
    void quicksort(index left, index right);

    /**
     * Partitions columnIdx between [@a left, @a right] after selecting the pivot in the middle.
     * @param left
     * @param right
     * @return The pivot.
     */
    index partition(index left, index right);

    /**
     * Binary search the sorted columnIdx vector between [@a left, @a right] for column @a j.
     * If @a j is not present, the index that is immediately left of the place where @a j would be
     * is returned. If
     * @param left
     * @param right
     * @param j
     * @return The position of column @a j in columnIdx or the element immediately to the left of the place where @a j
     * would be.
     */
    index binarySearchColumns(index left, index right, index j) const;

public:
    /** Default constructor */
    CSRMatrix();

    /**
     * Constructs the CSRMatrix with size @a dimension x @a dimension.
     * @param dimension Defines how many rows and columns this matrix has.
     * @param zero The zero element (default = 0.0).
     */
    CSRMatrix(count dimension, double zero = 0.0);

    /**
     * Constructs the CSRMatrix with size @a nRows x @a nCols.
     * @param nRows Number of rows.
     * @param nCols Number of columns.
     * @param zero The zero element (default = 0.0).
     */
    CSRMatrix(count nRows, count nCols, double zero = 0.0);

    /**
     * Constructs the @a dimension x @a dimension Matrix from the elements at position @a positions with values @values.
     * @param dimension Defines how many rows and columns this matrix has.
     * @param triplets The nonzero elements.
     * @param zero The zero element (default is 0.0).
     * @param isSorted True, if the triplets are sorted per row. Default is false.
     */
    CSRMatrix(count dimension, const std::vector<Triplet>& triplets, double zero = 0.0, bool isSorted = false);

    /**
     * Constructs the @a nRows x @a nCols Matrix from the elements at position @a positions with values @values.
     * @param nRows Defines how many rows this matrix has.
     * @param nCols Defines how many columns this matrix has.
     * @param triplets The nonzero elements.
     * @param zero The zero element (default is 0.0).
     * @param isSorted True, if the triplets are sorted per row. Default is false.
     */
    CSRMatrix(count nRows, count nCols, const std::vector<Triplet>& triplets, double zero = 0.0, bool isSorted = false);

    /**
     * Constructs the @a nRows x @a nCols Matrix from the elements stored in @a columnIdx and @a values. @a columnIdx and @a values store the colums and values by row.
     * @param nRows
     * @param nCols
     * @param columnIdx
     * @param values
     * @param zero The zero element (default is 0.0).
     * @param isSorted True if the column indices in @a columnIdx are sorted in every row.
     */
    CSRMatrix(count nRows, count nCols, const std::vector<std::vector<index>> &columnIdx, const std::vector<std::vector<double>> &values, double zero = 0.0, bool isSorted = false);

    /**
     * Constructs the @a nRows x @a nCols Matrix from the elements at position @a positions with values @values.
     * @param nRows Defines how many rows this matrix has.
     * @param nCols Defines how many columns this matrix has.
     * @param rowIdx The rowIdx vector of the CSR format.
     * @param columnIdx The columnIdx vector of the CSR format.
     * @param nonZeros The nonZero vector of the CSR format. Should be as long as the @a columnIdx vector.
     * @param zero The zero element (default is 0.0).
     * @param isSorted True, if the triplets are sorted per row. Default is false.
     */
    CSRMatrix(count nRows, count nCols, const std::vector<index>& rowIdx, const std::vector<index>& columnIdx, const std::vector<double>& nonZeros, double zero = 0.0, bool isSorted = false);

    /** Default copy constructor */
    CSRMatrix (const CSRMatrix &other) = default;

    /** Default move constructor */
    CSRMatrix (CSRMatrix &&other) = default;

    /** Default destructor */
    ~CSRMatrix() = default;

    /** Default move assignment operator */
    CSRMatrix& operator=(CSRMatrix &&other) = default;

    /** Default copy assignment operator */
    CSRMatrix& operator=(const CSRMatrix &other) = default;

    /**
     * Compares this matrix to @a other and returns true if the shape and zero element are the same as well as
     * all entries, otherwise returns false.
     * @param other
     */
    bool operator==(const CSRMatrix& other) const {
        bool equal = nRows == other.nRows && nCols == other.nCols && zero == other.zero;
        if (equal) {
            forNonZeroElementsInRowOrder([&](index i, index j, double value) {
                if (other(i,j) != value) {
                    equal = false;
                    return;
                }
            });
        }

        return equal;
    }

    /**
     * Compares this matrix to @a other and returns false if the shape and zero element are the same as well as
     * all entries, otherwise returns true.
     * @param other
     */
    bool operator!=(const CSRMatrix& other) const {
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
    count nnzInRow(index i) const;

    /**
     * @return Number of non-zeros in this matrix.
     */
    count nnz() const;

    /**
     * @return Value at matrix position (i,j).
     */
    double operator()(index i, index j) const;

    /**
     * Set the matrix at position (@a i, @a j) to @a value.
     * @note This operation can be linear in the number of non-zeros due to vector element movements
     */
    void setValue(index i, index j, double value);

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
    Vector row(index i) const;

    /**
     * @return Column @a j of this matrix as vector.
     */
    Vector column(index j) const;

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
    CSRMatrix operator*(double scalar) const;

    /**
     * Multiplies this matrix with a scalar specified in @a scalar.
     * @return Reference to this matrix.
     */
    CSRMatrix& operator*=(double scalar);

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
    CSRMatrix operator/(double divisor) const;

    /**
     * Divides this matrix by a divisor specified in @a divisor.
     * @return Reference to this matrix.
     */
    CSRMatrix& operator/=(double divisor);

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
     */
    CSRMatrix transpose() const;

    /**
     * Extracts a matrix with rows and columns specified by @a rowIndices and @a columnIndices from this matrix.
     * The order of rows and columns is equal to the order in @a rowIndices and @a columnIndices. It is also
     * possible to specify a row or column more than once to get duplicates.
     * @param rowIndices
     * @param columnIndices
     */
    CSRMatrix extract(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices) const;

    /**
     * Assign the contents of the matrix @a source to this matrix at rows and columns specified by @a rowIndices and
     * @a columnIndices. That is, entry (i,j) of @a source is assigned to entry (rowIndices[i], columnIndices[j]) of
     * this matrix. Note that the dimensions of @rowIndices and @a columnIndices must coincide with the number of rows
     * and columns of @a source.
     * @param rowIndices
     * @param columnIndices
     * @param source
     */
    void assign(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices, const CSRMatrix& source);

    /**
     * Applies the unary function @a unaryElementFunction to each value in the matrix. Note that it must hold that the
     * function applied to the zero element of this matrix returns the zero element.
     * @param unaryElementFunction
     */
    template<typename F>
    void apply(const F unaryElementFunction);

    /**
     * Compute the (weighted) adjacency matrix of the (weighted) Graph @a graph.
     * @param graph
     */
    static CSRMatrix adjacencyMatrix(const Graph& graph, double zero = 0.0);

    /**
     * Creates a diagonal matrix with dimension equal to the dimension of the Vector @a diagonalElements. The values on
     * the diagonal are the ones stored in @a diagonalElements (i.e. D(i,i) = diagonalElements[i]).
     * @param diagonalElements
     */
    static CSRMatrix diagonalMatrix(const Vector& diagonalElements, double zero = 0.0);

    /**
     * Returns the (weighted) incidence matrix of the (weighted) Graph @a graph.
     * @param graph
     */
    static CSRMatrix incidenceMatrix(const Graph& graph, double zero = 0.0);

    /**
     * Compute the (weighted) Laplacian of the (weighted) Graph @a graph.
     * @param graph
     */
    static CSRMatrix laplacianMatrix(const Graph& graph, double zero = 0.0);

    /**
     * Returns the (weighted) normalized Laplacian matrix of the (weighted) Graph @a graph
     * @param graph
     */
    static CSRMatrix normalizedLaplacianMatrix(const Graph& graph, double zero = 0.0);



    /**
     * Iterate over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
     */
    template<typename L> void forNonZeroElementsInRow(index row, L handle) const;

    /**
     * Iterate in parallel over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
     */
    template<typename L> void parallelForNonZeroElementsInRow(index row, L handle) const;

    /**
     * Iterate over all elements in row @a i in the matrix and call handle(index column, double value)
     */
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

template<typename L> inline CSRMatrix CSRMatrix::binaryOperator(const CSRMatrix &A, const CSRMatrix &B, L binaryOp) {
    assert(A.nRows == B.nRows && A.nCols == B.nCols);

    if (A.sorted() && B.sorted()) {
        std::vector<index> rowIdx(A.nRows+1);
        std::vector<std::vector<index>> columns(A.nRows);

        rowIdx[0] = 0;
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(A.nRows); ++i) {
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
        for (omp_index i = 0; i < static_cast<omp_index>(A.nRows); ++i) {
            for (index cIdx = rowIdx[i], j = 0; cIdx < rowIdx[i+1]; ++cIdx, ++j) {
                columnIdx[cIdx] = columns[i][j];
            }
            columns[i].clear();
            columns[i].resize(0);
            columns[i].shrink_to_fit();
        }

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(A.nRows); ++i) {
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

template<typename F>
void CSRMatrix::apply(F unaryElementFunction) {
#pragma omp parallel for
    for (omp_index k = 0; k < static_cast<omp_index>(nonZeros.size()); ++k) {
        nonZeros[k] = unaryElementFunction(nonZeros[k]);
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
    for (omp_index k = rowIdx[i]; k < static_cast<omp_index>(rowIdx[i+1]); ++k) {
        handle(columnIdx[k], nonZeros[k]);
    }
}

template<typename L>
inline void NetworKit::CSRMatrix::forElementsInRow(index i, L handle) const {
    Vector rowVector = row(i);
    index j = 0;
    rowVector.forElements([&](double val) {
        handle(j++, val);
    });
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
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
            handle(i, columnIdx[k], nonZeros[k]);
        }
    }
}

#endif // NETWORKIT_ALGEBRAIC_CSR_MATRIX_HPP_
