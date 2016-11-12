/*
 * DenseMatrixGTest.cpp
 *
 *  Created on: Nov 25, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "DenseMatrixGTest.h"

namespace NetworKit {

DenseMatrixGTest::DenseMatrixGTest() {
}

DenseMatrixGTest::~DenseMatrixGTest() {
}

TEST(DenseMatrixGTest, testMatrixDimension) {
	DenseMatrix mat(10, 10, std::vector<double>(100));

	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = DenseMatrix(5, 10, std::vector<double>(50));
	ASSERT_EQ(5u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = DenseMatrix(10, 5, std::vector<double>(50));
	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(5u, mat.numberOfColumns());
}


TEST(DenseMatrixGTest, testRowAndColumnAccess) {
	std::vector<double> values(100*100, 0.0);

	for (index i = 0; i < 100; ++i) {
		values[3 * 100 + i] = i;
	}

	values[10*100 + 10] = 42.123;

	DenseMatrix mat(100, 100, values);

	Vector v = mat.row(3);
	ASSERT_EQ(mat.numberOfColumns(), v.getDimension());

	for (index i = 0; i < 100; ++i) {
		EXPECT_EQ(i, v[i]);
	}

	v = mat.row(10);
	ASSERT_EQ(v.getDimension(), mat.numberOfColumns());
	ASSERT_TRUE(v.isTransposed());
	EXPECT_EQ(42.123, v[10]);

	v = mat.column(10);
	ASSERT_EQ(mat.numberOfRows(), v.getDimension());
	ASSERT_FALSE(v.isTransposed());

	EXPECT_EQ(10.0, v[3]);
	EXPECT_EQ(42.123, v[10]);

	mat.setValue(10, 10, 42);
	EXPECT_EQ(42, mat(10,10));


	// rectangular matrix
	// n x m (n < m)
	values = std::vector<double>(5*10, 0.0);

	values[4 * 10 + 9] = 11;
	values[0] = 42;

	mat = DenseMatrix(5, 10, values);
	v = mat.row(0);
	ASSERT_EQ(v.getDimension(), 10u);
	for (index i = 0; i < v.getDimension(); ++i) {
		if (i == 0) {
			EXPECT_EQ(42, v[i]);
		} else {
			EXPECT_EQ(0, v[i]);
		}
	}

	v = mat.column(9);
	ASSERT_EQ(v.getDimension(), 5u);
	for (index i = 0; i < v.getDimension(); ++i) {
		if (i == v.getDimension() - 1) {
			EXPECT_EQ(11, v[i]);
		} else {
			EXPECT_EQ(0, v[i]);
		}
	}

	// rectangular matrix
	// n x m (n > m)

	values = std::vector<double>(10*5, 0.0);

	values[9 * 5 + 4] = 11;
	values[0] = 42;

	mat = DenseMatrix(10, 5, values);
	v = mat.row(0);
	ASSERT_EQ(v.getDimension(), 5u);
	for (index i = 0; i < v.getDimension(); ++i) {
		if (i == 0) {
			EXPECT_EQ(42, v[i]);
		} else {
			EXPECT_EQ(0, v[i]);
		}
	}

	v = mat.column(4);
	ASSERT_EQ(v.getDimension(), 10u);
	for (index i = 0; i < v.getDimension(); ++i) {
		if (i == v.getDimension() - 1) {
			EXPECT_EQ(11, v[i]);
		} else {
			EXPECT_EQ(0, v[i]);
		}
	}
}

TEST(DenseMatrixGTest, testMatrixAddition) {
	std::vector<double> values1(100*100, 0.0);
	std::vector<double> values2(100*100, 0.0);

	for (index i = 0; i < 100; ++i) {
		values1[i * 100 + i] = 1;
		values2[i * 100 + i] = i;
	}

	values1[2 * 100 + 71] = 1.8;
	values2[42*100 + 43] = 3.14;

	DenseMatrix mat1(100, 100, values1);
	DenseMatrix mat2(100, 100, values2);

	DenseMatrix result = mat1 + mat2;
	ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

	EXPECT_EQ(0.0, result(10, 13));


	for (index i = 0; i < result.numberOfRows(); ++i) {
		EXPECT_EQ((i + 1), result(i, i));
	}
	EXPECT_EQ(1.8, result(2, 71));
	EXPECT_EQ(3.14, result(42, 43));

	EXPECT_EQ(0.0, result(3, 14));


//	// rectangular matrix
//	// n x m (n < m)
//	positions1 = {std::make_pair(0,0), std::make_pair(1,2)};
//	values1 = {1.0, 3.0};
//	mat1 = DenseMatrix(2, 5, positions1, values1);
//
//	positions2 = {std::make_pair(0,2), std::make_pair(1,2)};
//	values2 = {1.0, 1.0};
//	mat2 = DenseMatrix(2, 5, positions2, values2);
//
//	mat1.sort();
//	mat2.sort();
//
//	result = mat1 + mat2;
//
//	ASSERT_EQ(2u, result.numberOfRows());
//	ASSERT_EQ(5u, result.numberOfColumns());
//
//	EXPECT_EQ(1, result(0,0));
//	EXPECT_EQ(1, result(0,2));
//	EXPECT_EQ(4, result(1,2));
//
//	EXPECT_EQ(0, result(0,1));
//	EXPECT_EQ(0, result(1,4));
//
//	// rectangular matrix
//	// n x m (n > m)
//	positions1 = {std::make_pair(0,0), std::make_pair(2,1)};
//	values1 = {1.0, 3.0};
//	mat1 = DenseMatrix(5, 2, positions1, values1);
//
//	positions2 = {std::make_pair(2,0), std::make_pair(2,1)};
//	values2 = {1.0, 1.0};
//	mat2 = DenseMatrix(5, 2, positions2, values2);
//
//	mat1.sort();
//	mat2.sort();
//
//	result = mat1 + mat2;
//
//	ASSERT_EQ(5u, result.numberOfRows());
//	ASSERT_EQ(2u, result.numberOfColumns());
//
//	EXPECT_EQ(1, result(0,0));
//	EXPECT_EQ(1, result(2,0));
//	EXPECT_EQ(4, result(2,1));
//
//	EXPECT_EQ(0, result(0,1));
//	EXPECT_EQ(0, result(4,1));
}

//TEST(DenseMatrixGTest, testMatrixSubtraction) {
//	std::vector<std::pair<index, index> > positions1;
//	std::vector<std::pair<index, index> > positions2;
//	std::vector<double> values1;
//	std::vector<double> values2;
//
//	for (index i = 0; i < 100000; ++i) {
//		positions1.push_back(std::make_pair(i, i));
//		positions2.push_back(std::make_pair(i, i));
//		values1.push_back(1);
//		values2.push_back(i);
//	}
//
//	positions1.push_back(std::make_pair(2, 71));
//	values1.push_back(1.8);
//
//	positions2.push_back(std::make_pair(42, 43));
//	values2.push_back(3.14);
//
//	DenseMatrix mat1(100000, 100000, positions1, values1);
//	DenseMatrix mat2(100000, 100000, positions2, values2);
//
//	mat1.sort();
//	mat2.sort();
//
//	DenseMatrix result = mat2 - mat1;
//	ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
//	ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());
//
//	EXPECT_EQ(0.0, result(10, 13));
//
//	for (index i = 0; i < result.numberOfRows(); ++i) {
//		EXPECT_EQ(((int) i - 1), result(i, i));
//	}
//	EXPECT_EQ(-1.8, result(2, 71));
//	EXPECT_EQ(3.14, result(42, 43));
//
//	EXPECT_EQ(0.0, result(3, 14));
//
//	// rectangular matrix
//	// n x m (n < m)
//	positions1 = {std::make_pair(0,0), std::make_pair(1,2)};
//	values1 = {1.0, 3.0};
//	mat1 = DenseMatrix(2, 5, positions1, values1);
//
//
//	positions2 = {std::make_pair(0,2), std::make_pair(1,2)};
//	values2 = {1.0, 1.0};
//	mat2 = DenseMatrix(2, 5, positions2, values2);
//
//	mat1.sort();
//	mat2.sort();
//
//	result = mat1 - mat2;
//
//	ASSERT_EQ(2u, result.numberOfRows());
//	ASSERT_EQ(5u, result.numberOfColumns());
//
//	EXPECT_EQ(1, result(0,0));
//	EXPECT_EQ(-1, result(0,2));
//	EXPECT_EQ(2, result(1,2));
//
//	EXPECT_EQ(0, result(0,1));
//	EXPECT_EQ(0, result(1,4));
//
//	// rectangular matrix
//	// n x m (n > m)
//	positions1 = {std::make_pair(0,0), std::make_pair(2,1)};
//	values1 = {1.0, 3.0};
//	mat1 = DenseMatrix(5, 2, positions1, values1);
//
//	positions2 = {std::make_pair(2,0), std::make_pair(2,1)};
//	values2 = {1.0, 1.0};
//	mat2 = DenseMatrix(5, 2, positions2, values2);
//
//	mat1.sort();
//	mat2.sort();
//
//	result = mat1 - mat2;
//
//	ASSERT_EQ(5u, result.numberOfRows());
//	ASSERT_EQ(2u, result.numberOfColumns());
//
//	EXPECT_EQ(1, result(0,0));
//	EXPECT_EQ(-1, result(2,0));
//	EXPECT_EQ(2, result(2,1));
//
//	EXPECT_EQ(0, result(0,1));
//	EXPECT_EQ(0, result(4,1));
//}

//TEST(DenseMatrixGTest, testScalarMultiplication) {
//	std::vector<std::pair<index, index> > positions;
//	std::vector<double> values;
//
//	for (index i = 0; i < 10000; ++i) {
//		positions.push_back(std::make_pair(i, i));
//		values.push_back(i);
//	}
//
//	positions.push_back(std::make_pair(42, 43));
//	values.push_back(42.0);
//
//	DenseMatrix mat(10000, 10000, positions, values);
//	mat *= 2;
//	ASSERT_EQ(10000u, mat.numberOfRows());
//	ASSERT_EQ(10000u, mat.numberOfColumns());
//
//	for (index i = 0; i < 10000; ++i) {
//		EXPECT_EQ(i*2, mat(i, i));
//	}
//	EXPECT_EQ(84.0, mat(42, 43));
//	EXPECT_EQ(0.0, mat(55, 199));
//
//	mat *= 0.5;
//
//	for (index i = 0; i < 10000; ++i) {
//		EXPECT_EQ(i, mat(i, i));
//	}
//	EXPECT_EQ(42.0, mat(42, 43));
//	EXPECT_EQ(0.0, mat(55, 199));
//
//	// rectangular matrix
//	positions = {std::make_pair(0,0), std::make_pair(1,2)};
//	values = {1.0, 3.0};
//	mat = DenseMatrix(2, 5, positions, values);
//
//	mat *= 2;
//
//	EXPECT_EQ(2, mat(0,0));
//	EXPECT_EQ(6, mat(1,2));
//}

//TEST(DenseMatrixGTest, testMatrixDivisionOperator) {
//	std::vector<std::pair<index, index> > positions;
//	std::vector<double> values;
//
//	for (index i = 0; i < 10000; ++i) {
//		positions.push_back(std::make_pair(i, i));
//		values.push_back(i);
//	}
//
//	positions.push_back(std::make_pair(42, 43));
//	values.push_back(42.0);
//
//	DenseMatrix mat(10000, 10000, positions, values);
//	mat /= (1.0 / 2.0);
//	ASSERT_EQ(10000u, mat.numberOfRows());
//	ASSERT_EQ(10000u, mat.numberOfColumns());
//
//	for (index i = 0; i < 10000; ++i) {
//		EXPECT_EQ(i*2, mat(i, i));
//	}
//	EXPECT_EQ(84.0, mat(42, 43));
//	EXPECT_EQ(0.0, mat(55, 199));
//
//	mat /= 2;
//
//	for (index i = 0; i < 10000; ++i) {
//		EXPECT_EQ(i, mat(i, i));
//	}
//	EXPECT_EQ(42.0, mat(42, 43));
//	EXPECT_EQ(0.0, mat(55, 199));
//
//	// rectangular matrix
//	positions = {std::make_pair(0,0), std::make_pair(1,2)};
//	values = {1.0, 3.0};
//	mat = DenseMatrix(2, 5, positions, values);
//
//	mat /= 2;
//
//	EXPECT_EQ(0.5, mat(0,0));
//	EXPECT_EQ(1.5, mat(1,2));
//}
//
//TEST(DenseMatrixGTest, testMatrixVectorProduct) {
//	std::vector<std::pair<index, index> > mPositions;
//	std::vector<double> mValues;
//
//	for (index i = 0; i < 10000; ++i) {
//		mPositions.push_back(std::make_pair(i, i));
//		mValues.push_back(i);
//	}
//
//	mPositions.push_back(std::make_pair(42, 43));
//	mValues.push_back(42.0);
//
//	Vector vector(10000, 1.0);
//	vector[500] = 3.5;
//
//	DenseMatrix mat(10000, 10000, mPositions, mValues);
//
//	Vector result = mat * vector;
//	ASSERT_EQ(mat.numberOfRows(), result.getDimension());
//
//	for (index i = 0; i < 10000; ++i) {
//		if (i != 500 && i != 42 && i != 43) {
//			EXPECT_EQ(i, result[i]);
//		}
//	}
//
//	EXPECT_EQ(42.0, mat(42, 43));
//	EXPECT_EQ(84.0, result[42]);
//	EXPECT_EQ(1750.0, result[500]);
//
//
//	std::vector<std::pair<index, index> > positions;
//	positions.push_back(std::make_pair(0,0));
//	positions.push_back(std::make_pair(0,1));
//	positions.push_back(std::make_pair(0,2));
//	positions.push_back(std::make_pair(1,0));
//	positions.push_back(std::make_pair(1,1));
//	positions.push_back(std::make_pair(2,0));
//	positions.push_back(std::make_pair(2,2));
//	positions.push_back(std::make_pair(2,3));
//	positions.push_back(std::make_pair(3,2));
//	positions.push_back(std::make_pair(3,3));
//
//	std::vector<double> values = {1, 2, 3, 2, 2, 3, 3, -1, -1, 4};
//	DenseMatrix mat2(4, 4, positions, values);
//
//	Vector v({1,2,3,0});
//	Vector res = mat2 * v;
//	ASSERT_EQ(mat2.numberOfRows(), res.getDimension());
//
//	EXPECT_EQ(14, res[0]);
//	EXPECT_EQ(6, res[1]);
//	EXPECT_EQ(12, res[2]);
//	EXPECT_EQ(-3, res[3]);
//
//	// rectangular matrix
//	positions = {std::make_pair(0,0), std::make_pair(1,2)};
//	values = {1.0, 3.0};
//	mat = DenseMatrix(2, 5, positions, values);
//
//	v = {0,1,2,3,0};
//	res = mat * v;
//
//	ASSERT_EQ(2u, res.getDimension());
//	EXPECT_EQ(0, res[0]);
//	EXPECT_EQ(6, res[1]);
//}
//
//TEST(DenseMatrixGTest, testMatrixMultiplication) {
//	std::vector<std::pair<index, index> > positions;
//	std::vector<double> values = {1, 2, 3, 2, 2, 3, 3, -1, -1, 4};
//
//	positions.push_back(std::make_pair(0,0));
//	positions.push_back(std::make_pair(0,1));
//	positions.push_back(std::make_pair(0,2));
//	positions.push_back(std::make_pair(1,0));
//	positions.push_back(std::make_pair(1,1));
//	positions.push_back(std::make_pair(2,0));
//	positions.push_back(std::make_pair(2,2));
//	positions.push_back(std::make_pair(2,3));
//	positions.push_back(std::make_pair(3,2));
//	positions.push_back(std::make_pair(3,3));
//
//	//
//	//				 1  2  3  0
//	// 				 2  2  0  0
//	// mat1 = mat2 = 3  0  3 -1
//	//				 0  0 -1  4
//	//
//	DenseMatrix mat1(4, 4, positions, values);
//	ASSERT_EQ(4u, mat1.numberOfRows());
//	ASSERT_EQ(4u, mat1.numberOfColumns());
//
//	DenseMatrix mat2(4, 4, positions, values);
//	ASSERT_EQ(4u, mat2.numberOfRows());
//	ASSERT_EQ(4u, mat2.numberOfColumns());
//
//	//
//	//			14  6  12  -3
//	//			 6  8   6   0
//	// result = 12  6  19  -7
//	//			-3  0  -7  17
//	//
//	DenseMatrix result = mat1 * mat2;
//	ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
//	ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());
//	EXPECT_EQ(14u, result.nnz());
//
//	EXPECT_EQ(14, result(0,0));
//	EXPECT_EQ(6, result(0,1));
//	EXPECT_EQ(12, result(0,2));
//	EXPECT_EQ(-3, result(0,3));
//	EXPECT_EQ(6, result(1,0));
//	EXPECT_EQ(8, result(1,1));
//	EXPECT_EQ(6, result(1,2));
//	EXPECT_EQ(0, result(1,3));
//	EXPECT_EQ(12, result(2,0));
//	EXPECT_EQ(6, result(2,1));
//	EXPECT_EQ(19, result(2,2));
//	EXPECT_EQ(-7, result(2,3));
//	EXPECT_EQ(-3, result(3,0));
//	EXPECT_EQ(0, result(3,1));
//	EXPECT_EQ(-7, result(3,2));
//	EXPECT_EQ(17, result(3,3));
//
//
//	// rectangular matrices
//	positions = {std::make_pair(0,0), std::make_pair(0,3), std::make_pair(1,2), std::make_pair(2,1), std::make_pair(2,3)};
//	values = {1.0, 2.0, 1.0, 2.0, 4.0};
//	mat1 = DenseMatrix(3, 4, positions, values);
//
//	positions = {std::make_pair(0,0), std::make_pair(2,1), std::make_pair(3,0), std::make_pair(3,1)};
//	values = {1.0, 0.5, 42.0, 1.0};
//	mat2 = DenseMatrix(4, 2, positions, values);
//
//	result = mat1 * mat2;
//
//	EXPECT_EQ(85, result(0,0));
//	EXPECT_EQ(2, result(0,1));
//	EXPECT_EQ(0, result(1,0));
//	EXPECT_EQ(0.5, result(1,1));
//	EXPECT_EQ(168, result(2,0));
//	EXPECT_EQ(4, result(2,1));
//}
//
//TEST(DenseMatrixGTest, testBigMatrixMultiplication) {
//	METISGraphReader graphReader;
//	Graph G = graphReader.read("input/PGPgiantcompo.graph");
//
//	std::vector<std::pair<index,index>> positions;
//	std::vector<double> values;
//
//	G.forEdges([&](index i, index j, double value) {
//		positions.push_back(std::make_pair(i,j));
//		values.push_back(value);
//	});
//
//	DenseMatrix mat(G.upperNodeIdBound(), G.upperNodeIdBound(), positions, values);
//
//	DenseMatrix result = mat * mat;
//	ASSERT_EQ(mat.numberOfRows(), result.numberOfRows());
//	ASSERT_EQ(mat.numberOfColumns(), result.numberOfColumns());
//}
//
//TEST(DenseMatrixGTest, testTransposition) {
//	//
//	//	   1  2  3  1  1
//	// 	   0  2  0  0  0
//	// mat 4  0  3 -1  0
//	//	   0  0  0  4 -1
//	//
//	std::vector<std::pair<index, index> > positions;
//	std::vector<double> values = {1, 2, 3, 1, 1, 2, 4, 3, -1, 4, -1};
//
//	positions.push_back(std::make_pair(0,0));
//	positions.push_back(std::make_pair(0,1));
//	positions.push_back(std::make_pair(0,2));
//	positions.push_back(std::make_pair(0,3));
//	positions.push_back(std::make_pair(0,4));
//	positions.push_back(std::make_pair(1,1));
//	positions.push_back(std::make_pair(2,0));
//	positions.push_back(std::make_pair(2,2));
//	positions.push_back(std::make_pair(2,3));
//	positions.push_back(std::make_pair(3,3));
//	positions.push_back(std::make_pair(3,4));
//
//	DenseMatrix mat(4, 5, positions, values);
//	DenseMatrix matT = mat.transpose();
//
//	EXPECT_EQ(5u, matT.numberOfRows());
//	EXPECT_EQ(4u, matT.numberOfColumns());
//
//	mat.forNonZeroElementsInRowOrder([&](index i, index j, double value) {
//		EXPECT_EQ(value, matT(j,i));
//	});
//}
//
//TEST(DenseMatrixGTest, testMatrixTransposeMatrixMultiplication) {
//	std::vector<std::pair<index,index>> positions = {std::make_pair(0,0), std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,1), std::make_pair(2,0), std::make_pair(3,2)};
//	std::vector<double> values = {1.0, 2.0, 3.0, 2.0, 3.0, -1.0};
//	DenseMatrix A(4, 3, positions, values);
//
//	positions = {std::make_pair(0,0), std::make_pair(1,0), std::make_pair(2,1), std::make_pair(3,1), std::make_pair(3,2)};
//	values = {1.0, 3.0, -2.0, 5.0, -8.0};
//	DenseMatrix B(4, 3, positions, values);
//
//	DenseMatrix C = DenseMatrix::mTmMultiply(A, B);
//
//	EXPECT_EQ(1, C(0,0));
//	EXPECT_EQ(-6, C(0,1));
//	EXPECT_EQ(0, C(0,2));
//	EXPECT_EQ(8, C(1,0));
//	EXPECT_EQ(0, C(1,1));
//	EXPECT_EQ(0, C(1,2));
//	EXPECT_EQ(3, C(2,0));
//	EXPECT_EQ(-5, C(2,1));
//	EXPECT_EQ(8, C(2,2));
//}
//
//TEST(DenseMatrixGTest, testMatrixMatrixTransposeMultiplication) {
//	std::vector<std::pair<index,index>> positions = {std::make_pair(0,0), std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,1), std::make_pair(2,0), std::make_pair(3,2)};
//	std::vector<double> values = {1.0, 2.0, 3.0, 2.0, 3.0, -1.0};
//	DenseMatrix A(4, 3, positions, values);
//
//	positions = {std::make_pair(0,0), std::make_pair(1,0), std::make_pair(2,1), std::make_pair(3,1), std::make_pair(3,2)};
//	values = {1.0, 3.0, -2.0, 5.0, -8.0};
//	DenseMatrix B(4, 3, positions, values);
//
//	DenseMatrix C = DenseMatrix::mmTMultiply(A, B);
//
//	EXPECT_EQ(1, C(0,0));
//	EXPECT_EQ(3, C(0,1));
//	EXPECT_EQ(-4, C(0,2));
//	EXPECT_EQ(-14, C(0,3));
//	EXPECT_EQ(-4, C(1,2));
//	EXPECT_EQ(10, C(1,3));
//	EXPECT_EQ(3, C(2,0));
//	EXPECT_EQ(9, C(2,1));
//	EXPECT_EQ(8, C(3,3));
//	EXPECT_EQ(0, C(1,0));
//	EXPECT_EQ(0, C(1,1));
//	EXPECT_EQ(0, C(2,3));
//}
//
//TEST(DenseMatrixGTest, testMatrixTransposeVectorMultiplication) {
//	std::vector<std::pair<index,index>> positions = {std::make_pair(0,0), std::make_pair(2,1)};
//	std::vector<double> values = {1.0, 3.0};
//	DenseMatrix mat(5, 2, positions, values);
//
//	Vector v = {0,1,2,3,0};
//	Vector res = DenseMatrix::mTvMultiply(mat, v);
//
//	ASSERT_EQ(2u, res.getDimension());
//	EXPECT_EQ(0, res[0]);
//	EXPECT_EQ(6, res[1]);
//}
//
//TEST(DenseMatrixGTest, testMatrixDiagonal) {
//	//
//	//	   1  2  3  1  1
//	// 	   0  2  0  0  0
//	// mat 4  0  0 -1  0
//	//	   0  0  0  4 -1
//	//
//	std::vector<std::pair<index, index> > positions;
//	std::vector<double> values = {1, 2, 3, 1, 1, 2, 4, -1, 4, -1};
//
//	positions.push_back(std::make_pair(0,0));
//	positions.push_back(std::make_pair(0,1));
//	positions.push_back(std::make_pair(0,2));
//	positions.push_back(std::make_pair(0,3));
//	positions.push_back(std::make_pair(0,4));
//	positions.push_back(std::make_pair(1,1));
//	positions.push_back(std::make_pair(2,0));
//	positions.push_back(std::make_pair(2,3));
//	positions.push_back(std::make_pair(3,3));
//	positions.push_back(std::make_pair(3,4));
//
//	DenseMatrix mat(4, 5, positions, values);
//
//	Vector diag1 = mat.diagonal();
//	EXPECT_EQ(1, diag1[0]);
//	EXPECT_EQ(2, diag1[1]);
//	EXPECT_EQ(0, diag1[2]);
//	EXPECT_EQ(4, diag1[3]);
//
//	mat.sort();
//	Vector diag2 = mat.diagonal();
//	for (index i = 0; i < diag2.getDimension(); ++i) {
//		EXPECT_EQ(diag1[i], diag2[i]);
//	}
//}

TEST(DenseMatrixGTest, testLUDecomposition) {
	// 		  1  2  4
	// mat1 = 3  8  14
	// 		  2  6  13


	std::vector<double> values = {1,2,4,3,8,14,2,6,13};
	DenseMatrix mat(3, 3, values);
	DenseMatrix::LUDecomposition(mat);

	EXPECT_EQ(1, mat(0,0));
	EXPECT_EQ(2, mat(0,1));
	EXPECT_EQ(4, mat(0,2));
	EXPECT_EQ(3, mat(1,0));
	EXPECT_EQ(2, mat(1,1));
	EXPECT_EQ(2, mat(1,2));
	EXPECT_EQ(2, mat(2,0));
	EXPECT_EQ(1, mat(2,1));
	EXPECT_EQ(3, mat(2,2));

	Vector expected = {3, 4, -2};
	Vector b = {3, 13, 4};
	Vector luResult = DenseMatrix::LUSolve(mat, b);

	for (index i = 0; i < expected.getDimension(); ++i) {
		EXPECT_EQ(expected[i], luResult[i]);
	}
}

} /* namespace NetworKit */
