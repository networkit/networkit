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

namespace GraphBLAS {


// ****************************************************
// 						Operations
// ****************************************************

template<class SEMIRING = ArithmeticSemiring, class MATRIX>
MATRIX MxM(const MATRIX& A, const MATRIX& B) {
	assert(A.numberOfColumns() == B.numberOfRows());
	assert(A.getZero() == SEMIRING::zero() && B.getZero() == SEMIRING::zero());

	std::vector<NetworKit::Triplet> triplets;
	NetworKit::SparseAccumulator spa(A.numberOfRows());
	for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
		A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
			B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
				spa.scatter(SEMIRING::mult(w1,w2), j, [](const double a, const double b) {return SEMIRING::add(a,b);});
			});
		});

		spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
			triplets.push_back({i,j,value});
		});

		spa.increaseRow();
	}

	return MATRIX(A.numberOfRows(), B.numberOfColumns(), triplets, A.getZero());
}

//template<class SEMIRING = ArithmeticSemiring, class MATRIX>
//void MxM(const MATRIX& A, const MATRIX& B, MATRIX& C) {
//	assert(A.numberOfColumns() == B.numberOfRows() && A.numberOfRows() == C.numberOfRows() && B.numberOfColumns() = C.numberOfColumns());
//	assert(A.getZero() == SEMIRING::zero() && B.getZero() == SEMIRING::getZero() && C.getZero() == SEMIRING::getZero());
//
//	NetworKit::SparseAccumulator spa(A.numberOfRows());
//	for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
//		A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
//			B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
//				spa.scatter(SEMIRING::mult(w1,w2), j, [](const double a, const double b) {return SEMIRING::add(a,b);});
//			});
//		});
//
//		spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
//			C.setValue(i,j, C(i,j) + value);
//		});
//
//		spa.increaseRow();
//	}
//}

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

template<class SEMIRING, class MATRIX, typename L>
MATRIX eWiseBinOp(const MATRIX& A, const MATRIX& B, L binOp) {
	assert(A.numberOfRows() == B.numberOfRows() && A.numberOfColumns() == B.numberOfColumns());
	assert(A.getZero() == B.getZero());

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


template<class SEMIRING = ArithmeticSemiring, class MATRIX>
MATRIX eWiseAdd(const MATRIX& A, const MATRIX& B) {
	return eWiseBinOp<SEMIRING, MATRIX>(A, B, [](const double a, const double b) {return SEMIRING::add(a,b);});
}

template<class SEMIRING = ArithmeticSemiring, class MATRIX>
MATRIX eWiseMult(const MATRIX& A, const MATRIX& B) {
	return eWiseBinOp<SEMIRING, MATRIX>(A, B, [](const double a, const double b) {return SEMIRING::mult(a,b);});
}













}



#endif /* NETWORKIT_CPP_ALGEBRAIC_GRAPHBLAS_H_ */
