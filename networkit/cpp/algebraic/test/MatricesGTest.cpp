/*
 * MatricesGTest.cpp
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MatricesGTest.h"

#include "../CSRMatrix.h"
#include "../DenseMatrix.h"
#include "../DynamicMatrix.h"

namespace NetworKit {

TEST_F(MatricesGTest, testDimension) {
	testDimension<DynamicMatrix>();
	testDimension<CSRMatrix>();
	testDimension<DenseMatrix>();
}

TEST_F(MatricesGTest, testNNZInRow) {
	testNNZInRow<DynamicMatrix>();
	testNNZInRow<CSRMatrix>();
	testNNZInRow<DenseMatrix>();
}

TEST_F(MatricesGTest, testRowAndColumnAccess) {
	testRowAndColumnAccess<DynamicMatrix>();
	testRowAndColumnAccess<CSRMatrix>();
	testRowAndColumnAccess<DenseMatrix>();
}

TEST_F(MatricesGTest, testDiagonalVector) {
	testDiagonalVector<DynamicMatrix>();
	testDiagonalVector<CSRMatrix>();
	testDiagonalVector<DenseMatrix>();
}

TEST_F(MatricesGTest, testTranspose) {
	testTranspose<DynamicMatrix>();
	testTranspose<CSRMatrix>();
	testTranspose<DenseMatrix>();
}

TEST_F(MatricesGTest, testExtract) {
	testExtract<DynamicMatrix>();
	testExtract<CSRMatrix>();
}

TEST_F(MatricesGTest, testAssign) {
	testAssign<DynamicMatrix>();
	testAssign<CSRMatrix>();
	testAssign<DenseMatrix>();
}

TEST_F(MatricesGTest, testApply) {
	testApply<DynamicMatrix>();
	testApply<CSRMatrix>();
}

TEST_F(MatricesGTest, testMatrixAddition) {
	testMatrixAddition<DynamicMatrix>();
	testMatrixAddition<CSRMatrix>();
	testMatrixAddition<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixSubtraction) {
	testMatrixSubtraction<DynamicMatrix>();
	testMatrixSubtraction<CSRMatrix>();
	testMatrixSubtraction<DenseMatrix>();
}

TEST_F(MatricesGTest, testScalarMultiplication) {
	testScalarMultiplication<DynamicMatrix>();
	testScalarMultiplication<CSRMatrix>();
	testScalarMultiplication<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixDivisionOperator) {
	testMatrixDivisionOperator<DynamicMatrix>();
	testMatrixDivisionOperator<CSRMatrix>();
	testMatrixDivisionOperator<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixVectorProduct) {
	testMatrixVectorProduct<DynamicMatrix>();
	testMatrixVectorProduct<CSRMatrix>();
	testMatrixVectorProduct<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixMultiplication) {
	testMatrixMultiplication<DynamicMatrix>();
	testMatrixMultiplication<CSRMatrix>();
	testMatrixMultiplication<DenseMatrix>();
}

TEST_F(MatricesGTest, testBigMatrixMultiplcation) {
	testBigMatrixMultiplication<DynamicMatrix>();
	testBigMatrixMultiplication<CSRMatrix>();
}

TEST_F(MatricesGTest, testAdjacencyMatrixOfGraph) {
	testAdjacencyMatrix<DynamicMatrix>();
	testAdjacencyMatrix<CSRMatrix>();
}

TEST_F(MatricesGTest, testDiagonalMatrix) {
	testDiagonalMatrix<DynamicMatrix>();
	testDiagonalMatrix<CSRMatrix>();
}

TEST_F(MatricesGTest, testIncidenceMatrix) {
	testIncidenceMatrix<DynamicMatrix>();
	testIncidenceMatrix<CSRMatrix>();
}

TEST_F(MatricesGTest, testLaplacianMatrixOfGraph) {
	testLaplacianOfGraph<DynamicMatrix>();
	testLaplacianOfGraph<CSRMatrix>();
}

} /* namespace NetworKit */
