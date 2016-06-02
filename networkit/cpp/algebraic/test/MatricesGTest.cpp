/*
 * MatricesGTest.cpp
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MatricesGTest.h"

#include "../Matrix.h"
#include "../CSRMatrix.h"
#include "../DenseMatrix.h"

namespace NetworKit {

TEST_F(MatricesGTest, testDimension) {
	testDimension<Matrix>();
	testDimension<CSRMatrix>();
	testDimension<DenseMatrix>();
}

TEST_F(MatricesGTest, testNNZInRow) {
	testNNZInRow<Matrix>();
	testNNZInRow<CSRMatrix>();
	testNNZInRow<DenseMatrix>();
}

TEST_F(MatricesGTest, testRowAndColumnAccess) {
	testRowAndColumnAccess<Matrix>();
	testRowAndColumnAccess<CSRMatrix>();
	testRowAndColumnAccess<DenseMatrix>();
}

TEST_F(MatricesGTest, testDiagonalVector) {
	testDiagonalVector<Matrix>();
	testDiagonalVector<CSRMatrix>();
	testDiagonalVector<DenseMatrix>();
}

TEST_F(MatricesGTest, testTranspose) {
	testTranspose<Matrix>();
	testTranspose<CSRMatrix>();
	testTranspose<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixAddition) {
	testMatrixAddition<Matrix>();
	testMatrixAddition<CSRMatrix>();
	testMatrixAddition<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixSubtraction) {
	testMatrixSubtraction<Matrix>();
	testMatrixSubtraction<CSRMatrix>();
	testMatrixSubtraction<DenseMatrix>();
}

TEST_F(MatricesGTest, testScalarMultiplication) {
	testScalarMultiplication<Matrix>();
	testScalarMultiplication<CSRMatrix>();
	testScalarMultiplication<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixDivisionOperator) {
	testMatrixDivisionOperator<Matrix>();
	testMatrixDivisionOperator<CSRMatrix>();
	testMatrixDivisionOperator<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixVectorProduct) {
	testMatrixVectorProduct<Matrix>();
	testMatrixVectorProduct<CSRMatrix>();
	testMatrixVectorProduct<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixMultiplication) {
	testMatrixMultiplication<Matrix>();
	testMatrixMultiplication<CSRMatrix>();
	testMatrixMultiplication<DenseMatrix>();
}

TEST_F(MatricesGTest, testBigMatrixMultiplcation) {
	testBigMatrixMultiplication<Matrix>();
	testBigMatrixMultiplication<CSRMatrix>();
}

TEST_F(MatricesGTest, testLaplacianMatrixOfGraph) {
	testLaplacianOfGraph<Matrix>();
	testLaplacianOfGraph<CSRMatrix>();
}

} /* namespace NetworKit */
