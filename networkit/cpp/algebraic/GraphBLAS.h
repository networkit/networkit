/*
 * GraphBLAS.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_GRAPHBLAS_H_
#define NETWORKIT_CPP_ALGEBRAIC_GRAPHBLAS_H_

#include <limits>
#include "Semirings.h"
#include "SparseAccumulator.h"
#include "AlgebraicGlobals.h"
#include "Vector.h"

/**
 * @ingroup algebraic
 */
namespace GraphBLAS {

// ****************************************************
// 						Operations
// ****************************************************

/**
 * Computes binOp(A(i,j), B(i,j)) for all i,j element-wise. Note that the dimensions of
 * @a A and @a B must coincide and that the zero must be the same.
 * @param A
 * @param B
 * @param binOp
 * @return The resulting matrix.
 */
template<class SEMIRING, class MATRIX, typename L>
MATRIX eWiseBinOp(const MATRIX& A, const MATRIX& B, L binOp) {
	assert(A.numberOfRows() == B.numberOfRows() && A.numberOfColumns() == B.numberOfColumns());
	assert(A.getZero() == B.getZero() && A.getZero() == SEMIRING::zero());

	std::vector<int64_t> columnPointer(A.numberOfColumns(), -1);
	std::vector<double> Arow(A.numberOfColumns(), SEMIRING::zero());
	std::vector<double> Brow(A.numberOfColumns(), SEMIRING::zero());

	std::vector<NetworKit::Triplet> triplets;

	for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
		NetworKit::index listHead = 0;
		NetworKit::count nnz = 0;

		// search for nonZeros in matrix A
		A.forNonZeroElementsInRow(i, [&](NetworKit::index j, double value) {
			Arow[j] = value;

			columnPointer[j] = listHead;
			listHead = j;
			nnz++;
		});

		// search for nonZeros in matrix B
		B.forNonZeroElementsInRow(i, [&](NetworKit::index j, double value) {
			Brow[j] = value;

			if (columnPointer[j] == -1) { // matrix A does not have a nonZero entry in column j
				columnPointer[j] = listHead;
				listHead = j;
				nnz++;
			}
		});

		// apply operator on the found nonZeros in A and B
		for (NetworKit::count k = 0; k < nnz; ++k) {
			double value = binOp(Arow[listHead], Brow[listHead]);
			if (value != SEMIRING::zero()) {
				triplets.push_back({i,listHead,value});
			}

			NetworKit::index temp = listHead;
			listHead = columnPointer[listHead];

			// reset for next row
			columnPointer[temp] = -1;
			Arow[temp] = SEMIRING::zero();
			Brow[temp] = SEMIRING::zero();
		}

		nnz = 0;
	}

	return MATRIX(A.numberOfRows(), A.numberOfColumns(), triplets, A.getZero());
}

/**
 * Computes the matrix-matrix multiplication of @a A and @a B. Note that
 * A.numberOfColumns() must be equal to B.numberOfRows() and the zero elements
 * must be the same. The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @return The result of the multiplication A * B.
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
MATRIX MxM(const MATRIX& A, const MATRIX& B) {
	assert(A.numberOfColumns() == B.numberOfRows());
	assert(A.getZero() == SEMIRING::zero() && B.getZero() == SEMIRING::zero());

	std::vector<NetworKit::Triplet> triplets;
	NetworKit::SparseAccumulator spa(B.numberOfRows());
	for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
		A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
			B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
				spa.scatter(SEMIRING::mult(w1,w2), j, *SEMIRING::add);
			});
		});

		spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
			triplets.push_back({i,j,value});
		});

		spa.increaseRow();
	}

	return MATRIX(A.numberOfRows(), B.numberOfColumns(), triplets, A.getZero());
}

/**
 * Computes the matrix-matrix multiplication of @a A and @a B and adds it to @a C where
 * the add operation is that of the specified Semiring (i.e. C(i,j) = SEMIRING::add(C(i,j), (A*B)(i,j))).
 * The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @param C
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
void MxM(const MATRIX& A, const MATRIX& B, MATRIX& C) {
	assert(A.numberOfColumns() == B.numberOfRows() && A.numberOfRows() == C.numberOfRows() && B.numberOfColumns() == C.numberOfColumns());
	assert(A.getZero() == SEMIRING::zero() && B.getZero() == SEMIRING::zero() && C.getZero() == SEMIRING::zero());

	std::vector<NetworKit::Triplet> triplets;
	NetworKit::SparseAccumulator spa(B.numberOfRows());
	for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
		A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
			B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
				spa.scatter(SEMIRING::mult(w1,w2), j, *SEMIRING::add);
			});
		});

		spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
			triplets.push_back({i,j,value});
		});

		spa.increaseRow();
	}

	MATRIX temp(A.numberOfRows(), B.numberOfRows(), triplets, A.getZero());
	C = eWiseBinOp<SEMIRING, MATRIX>(C, temp, *SEMIRING::add);
}

/**
 * Computes the matrix-matrix multiplication of @a A and @a B and adds it to @a C where
 * the add operation is specified by the binary function @a accum (i.e. C(i,j) = accum(C(i,j), (A*B)(i,j))).
 * The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @param C
 * @param accum
 */
template<class SEMIRING = ArithmeticSemiring, typename F, class MATRIX>
void MxM(const MATRIX& A, const MATRIX& B, MATRIX& C, F accum) {
	assert(A.numberOfColumns() == B.numberOfRows() && A.numberOfRows() == C.numberOfRows() && B.numberOfColumns() == C.numberOfColumns());
	assert(A.getZero() == SEMIRING::zero() && B.getZero() == SEMIRING::zero() && C.getZero() == SEMIRING::zero());

	std::vector<NetworKit::Triplet> triplets;
	NetworKit::SparseAccumulator spa(B.numberOfRows());
	for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
		A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
			B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
				spa.scatter(SEMIRING::mult(w1,w2), j, *SEMIRING::add);
			});
		});

		spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
			triplets.push_back({i,j,value});
		});

		spa.increaseRow();
	}

	MATRIX temp(A.numberOfRows(), B.numberOfRows(), triplets, A.getZero());
	C = eWiseBinOp<SEMIRING, MATRIX>(C, temp, accum);
}

/**
 * Computes the matrix-vector product of matrix @a A and Vector @a v. The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param v
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
NetworKit::Vector MxV(const MATRIX& A, const NetworKit::Vector& v) {
	assert(!v.isTransposed());
	assert(A.numberOfColumns() == v.getDimension());
	assert(A.getZero() == SEMIRING::zero());
	NetworKit::Vector result(A.numberOfRows(), A.getZero());

	A.parallelForNonZeroElementsInRowOrder([&](NetworKit::index i, NetworKit::index j, double value) {
		result[i] = SEMIRING::add(result[i], SEMIRING::mult(value, v[j]));
	});

	return result;
}

/**
 * Computes the matrix-vector product of matrix @a A and Vector @a v and adds it to @a c where the add operation
 * is that of the specified Semiring (i.e. c[i] = SEMIRING::add(c[i], (A*v)[i]). The default Semiring is the
 * ArithmeticSemiring.
 * @param A
 * @param v
 * @param c
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
void MxV(const MATRIX& A, const NetworKit::Vector& v, NetworKit::Vector& c) {
	assert(!v.isTransposed());
	assert(A.numberOfColumns() == v.getDimension());
	assert(A.getZero() == SEMIRING::zero());

	A.parallelForNonZeroElementsInRowOrder([&](NetworKit::index i, NetworKit::index j, double value) {
		c[i] = SEMIRING::add(c[i], SEMIRING::mult(value, v[j]));
	});
}

/**
 * Computes the matrix-vector product of matrix @a A and Vector @a v and adds it to @a c where the add operation
 * is that of the specified binary function @a accum (i.e. c[i] = accum(c[i], (A*v)[i]). The default Semiring is the
 * ArithmeticSemiring.
 * @param A
 * @param v
 * @param c
 */
template<class SEMIRING = ArithmeticSemiring, typename F, class MATRIX>
void MxV(const MATRIX& A, const NetworKit::Vector& v, NetworKit::Vector& c, F accum) {
	assert(!v.isTransposed());
	assert(A.numberOfColumns() == v.getDimension());
	assert(A.getZero() == SEMIRING::zero());

	A.parallelForNonZeroElementsInRowOrder([&](NetworKit::index i, NetworKit::index j, double value) {
		c[i] = accum(c[i], SEMIRING::mult(value, v[j]));
	});
}

/**
 * Computes SEMIRING::add(A(i,j), B(i,j)) for all i,j element-wise and returns the resulting matrix. The default
 * Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
MATRIX eWiseAdd(const MATRIX& A, const MATRIX& B) {
	return eWiseBinOp<SEMIRING, MATRIX>(A, B, [](const double a, const double b) {return SEMIRING::add(a,b);});
}

/**
 * Computes SEMIRING::mult(A(i,j), B(i,j)) for all i,j element-wise and returns the resulting matrix. The default
 * Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @return
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
MATRIX eWiseMult(const MATRIX& A, const MATRIX& B) {
	return eWiseBinOp<SEMIRING, MATRIX>(A, B, [](const double a, const double b) {return SEMIRING::mult(a,b);});
}

/**
 * Computes the row-reduction of the @a matrix and returns the result as a vector. That is, the elements of each row
 * are summed up to form the respective entry in the result vector. The add operator is that of the specified
 * Semiring. The default Semiring is the ArithmeticSemiring.
 * @param matrix
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
NetworKit::Vector rowReduce(const MATRIX& matrix) {
	assert(matrix.getZero() == SEMIRING::zero());
	NetworKit::Vector rowReduction(matrix.numberOfRows(), 0.0);

#pragma omp parallel for
	for (NetworKit::index i = 0; i < matrix.numberOfRows(); ++i) {
		matrix.forNonZeroElementsInRow(i, [&](NetworKit::index j, double value) {
			rowReduction[i] = SEMIRING::add(rowReduction[i], value);
		});
	}

	return rowReduction;
}

/**
 * Computes the column-reduction of the @a matrix and returns the result as a Vector. That is, the elements of each
 * column are summed up to form the respective entry in the result Vector. The add operator is that of the specified
 * Semiring. The default Semiring is the ArithmeticSemiring.
 * @param matrix
 */
template<class SEMIRING = ArithmeticSemiring, class MATRIX>
NetworKit::Vector columnReduce(const MATRIX& matrix) {
	assert(matrix.getZero() == SEMIRING::zero());
	NetworKit::Vector columnReduction(matrix.numberOfColumns(), 0.0);

	matrix.forNonZeroElementsInRowOrder([&](NetworKit::index i, NetworKit::index j, double value) {
		columnReduction[j] = SEMIRING::add(columnReduction[j], value);
	});

	return columnReduction;
}

}



#endif /* NETWORKIT_CPP_ALGEBRAIC_GRAPHBLAS_H_ */
