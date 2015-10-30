/*
 * CSRMatrixGTest.cpp
 *
 *  Created on: May 13, 2015
 *      Author: Michael
 */

#include "CSRMatrixGTest.h"

namespace NetworKit {

CSRMatrixGTest::CSRMatrixGTest() {
}

CSRMatrixGTest::~CSRMatrixGTest() {
}

TEST(CSRMatrixGTest, testMatrixDimension) {
	CSRMatrix mat(10, 10, std::vector<std::pair<index,index>>(), std::vector<double>());;

	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = CSRMatrix (5, 10, std::vector<std::pair<index,index>>(), std::vector<double>());;
	ASSERT_EQ(5u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = CSRMatrix (10, 5, std::vector<std::pair<index,index>>(), std::vector<double>());;
	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(5u, mat.numberOfColumns());
}

TEST(CSRMatrixGTest, testNNZInRow) {
	std::vector<std::pair<index,index>> positions = {std::make_pair(0,0), std::make_pair(0,2), std::make_pair(1,0), std::make_pair(3,2)};
	std::vector<double> values = {1.0, 2.0, 4.0, 2.0};

	CSRMatrix mat(4, 3, positions, values);
	EXPECT_EQ(2u, mat.nnzInRow(0));
	EXPECT_EQ(1u, mat.nnzInRow(1));
	EXPECT_EQ(0u, mat.nnzInRow(2));
	EXPECT_EQ(1u, mat.nnzInRow(3));
}

TEST(CSRMatrixGTest, testRowAndColumnAccess) {
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values;

	for (index i = 0; i < 1000; ++i) {
		positions.push_back(std::make_pair(3, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(10, 10));
	values.push_back(42.123);

	CSRMatrix mat(1000, 1000, positions, values);

	Vector v = mat.row(3);
	ASSERT_EQ(mat.numberOfColumns(), v.getDimension());

	for (index i = 0; i < 1000; ++i) {
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


	// rectangular matrix
	// n x m (n < m)
	positions.clear();
	values.clear();

	positions.push_back(std::make_pair(4,9));
	values.push_back(11);

	positions.push_back(std::make_pair(0,0));
	values.push_back(42);

	mat = CSRMatrix(5, 10, positions, values);
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

	positions.clear();
	values.clear();

	positions.push_back(std::make_pair(9,4));
	values.push_back(11);

	positions.push_back(std::make_pair(0,0));
	values.push_back(42);

	mat = CSRMatrix(10, 5, positions, values);
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

TEST(CSRMatrixGTest, testMatrixAddition) {
	std::vector<std::pair<index, index> > positions1;
	std::vector<std::pair<index, index> > positions2;
	std::vector<double> values1;
	std::vector<double> values2;

	for (index i = 0; i < 100000; ++i) {
		positions1.push_back(std::make_pair(i, i));
		positions2.push_back(std::make_pair(i, i));
		values1.push_back(1);
		values2.push_back(i);
	}

	positions1.push_back(std::make_pair(2, 71));
	values1.push_back(1.8);

	positions2.push_back(std::make_pair(42, 43));
	values2.push_back(3.14);

	CSRMatrix mat1(100000, 100000, positions1, values1);
	CSRMatrix mat2(100000, 100000, positions2, values2);

	CSRMatrix result = mat1 + mat2;
	ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

	EXPECT_EQ(0.0, result(10, 13));


	for (index i = 0; i < result.numberOfRows(); ++i) {
		EXPECT_EQ((i + 1), result(i, i));
	}
	EXPECT_EQ(1.8, result(2, 71));
	EXPECT_EQ(3.14, result(42, 43));

	EXPECT_EQ(0.0, result(3, 14));


	// rectangular matrix
	// n x m (n < m)
	positions1 = {std::make_pair(0,0), std::make_pair(1,2)};
	values1 = {1.0, 3.0};
	mat1 = CSRMatrix(2, 5, positions1, values1);

	positions2 = {std::make_pair(0,2), std::make_pair(1,2)};
	values2 = {1.0, 1.0};
	mat2 = CSRMatrix(2, 5, positions2, values2);

	result = mat1 + mat2;

	ASSERT_EQ(2u, result.numberOfRows());
	ASSERT_EQ(5u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(1, result(0,2));
	EXPECT_EQ(4, result(1,2));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(1,4));

	// rectangular matrix
	// n x m (n > m)
	positions1 = {std::make_pair(0,0), std::make_pair(2,1)};
	values1 = {1.0, 3.0};
	mat1 = CSRMatrix(5, 2, positions1, values1);

	positions2 = {std::make_pair(2,0), std::make_pair(2,1)};
	values2 = {1.0, 1.0};
	mat2 = CSRMatrix(5, 2, positions2, values2);

	result = mat1 + mat2;

	ASSERT_EQ(5u, result.numberOfRows());
	ASSERT_EQ(2u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(1, result(2,0));
	EXPECT_EQ(4, result(2,1));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(4,1));
}

TEST(CSRMatrixGTest, testMatrixSubtraction) {
	std::vector<std::pair<index, index> > positions1;
	std::vector<std::pair<index, index> > positions2;
	std::vector<double> values1;
	std::vector<double> values2;

	for (index i = 0; i < 100000; ++i) {
		positions1.push_back(std::make_pair(i, i));
		positions2.push_back(std::make_pair(i, i));
		values1.push_back(1);
		values2.push_back(i);
	}

	positions1.push_back(std::make_pair(2, 71));
	values1.push_back(1.8);

	positions2.push_back(std::make_pair(42, 43));
	values2.push_back(3.14);

	CSRMatrix mat1(100000, 100000, positions1, values1);
	CSRMatrix mat2(100000, 100000, positions2, values2);

	CSRMatrix result = mat2 - mat1;
	ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

	EXPECT_EQ(0.0, result(10, 13));

	for (index i = 0; i < result.numberOfRows(); ++i) {
		EXPECT_EQ(((int) i - 1), result(i, i));
	}
	EXPECT_EQ(-1.8, result(2, 71));
	EXPECT_EQ(3.14, result(42, 43));

	EXPECT_EQ(0.0, result(3, 14));

	// rectangular matrix
	// n x m (n < m)
	positions1 = {std::make_pair(0,0), std::make_pair(1,2)};
	values1 = {1.0, 3.0};
	mat1 = CSRMatrix(2, 5, positions1, values1);


	positions2 = {std::make_pair(0,2), std::make_pair(1,2)};
	values2 = {1.0, 1.0};
	mat2 = CSRMatrix(2, 5, positions2, values2);

	result = mat1 - mat2;

	ASSERT_EQ(2u, result.numberOfRows());
	ASSERT_EQ(5u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(-1, result(0,2));
	EXPECT_EQ(2, result(1,2));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(1,4));

	// rectangular matrix
	// n x m (n > m)
	positions1 = {std::make_pair(0,0), std::make_pair(2,1)};
	values1 = {1.0, 3.0};
	mat1 = CSRMatrix(5, 2, positions1, values1);

	positions2 = {std::make_pair(2,0), std::make_pair(2,1)};
	values2 = {1.0, 1.0};
	mat2 = CSRMatrix(5, 2, positions2, values2);

	result = mat1 - mat2;

	ASSERT_EQ(5u, result.numberOfRows());
	ASSERT_EQ(2u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(-1, result(2,0));
	EXPECT_EQ(2, result(2,1));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(4,1));
}

TEST(CSRMatrixGTest, testScalarMultiplication) {
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values;

	for (index i = 0; i < 10000; ++i) {
		positions.push_back(std::make_pair(i, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(42, 43));
	values.push_back(42.0);

	CSRMatrix mat(10000, 10000, positions, values);
	mat *= 2;
	ASSERT_EQ(10000u, mat.numberOfRows());
	ASSERT_EQ(10000u, mat.numberOfColumns());

	for (index i = 0; i < 10000; ++i) {
		EXPECT_EQ(i*2, mat(i, i));
	}
	EXPECT_EQ(84.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 199));

	mat *= 0.5;

	for (index i = 0; i < 10000; ++i) {
		EXPECT_EQ(i, mat(i, i));
	}
	EXPECT_EQ(42.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 199));

	// rectangular matrix
	positions = {std::make_pair(0,0), std::make_pair(1,2)};
	values = {1.0, 3.0};
	mat = CSRMatrix(2, 5, positions, values);

	mat *= 2;

	EXPECT_EQ(2, mat(0,0));
	EXPECT_EQ(6, mat(1,2));
}

TEST(CSRMatrixGTest, testMatrixDivisionOperator) {
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values;

	for (index i = 0; i < 10000; ++i) {
		positions.push_back(std::make_pair(i, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(42, 43));
	values.push_back(42.0);

	CSRMatrix mat(10000, 10000, positions, values);
	mat /= (1.0 / 2.0);
	ASSERT_EQ(10000u, mat.numberOfRows());
	ASSERT_EQ(10000u, mat.numberOfColumns());

	for (index i = 0; i < 10000; ++i) {
		EXPECT_EQ(i*2, mat(i, i));
	}
	EXPECT_EQ(84.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 199));

	mat /= 2;

	for (index i = 0; i < 10000; ++i) {
		EXPECT_EQ(i, mat(i, i));
	}
	EXPECT_EQ(42.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 199));

	// rectangular matrix
	positions = {std::make_pair(0,0), std::make_pair(1,2)};
	values = {1.0, 3.0};
	mat = CSRMatrix(2, 5, positions, values);

	mat /= 2;

	EXPECT_EQ(0.5, mat(0,0));
	EXPECT_EQ(1.5, mat(1,2));
}

TEST(CSRMatrixGTest, testMatrixVectorProduct) {
	std::vector<std::pair<index, index> > mPositions;
	std::vector<double> mValues;

	for (index i = 0; i < 10000; ++i) {
		mPositions.push_back(std::make_pair(i, i));
		mValues.push_back(i);
	}

	mPositions.push_back(std::make_pair(42, 43));
	mValues.push_back(42.0);

	Vector vector(10000, 1.0);
	vector[500] = 3.5;

	CSRMatrix mat(10000, 10000, mPositions, mValues);

	Vector result = mat * vector;
	ASSERT_EQ(mat.numberOfRows(), result.getDimension());

	for (index i = 0; i < 10000; ++i) {
		if (i != 500 && i != 42 && i != 43) {
			EXPECT_EQ(i, result[i]);
		}
	}

	EXPECT_EQ(42.0, mat(42, 43));
	EXPECT_EQ(84.0, result[42]);
	EXPECT_EQ(1750.0, result[500]);


	std::vector<std::pair<index, index> > positions;
	positions.push_back(std::make_pair(0,0));
	positions.push_back(std::make_pair(0,1));
	positions.push_back(std::make_pair(0,2));
	positions.push_back(std::make_pair(1,0));
	positions.push_back(std::make_pair(1,1));
	positions.push_back(std::make_pair(2,0));
	positions.push_back(std::make_pair(2,2));
	positions.push_back(std::make_pair(2,3));
	positions.push_back(std::make_pair(3,2));
	positions.push_back(std::make_pair(3,3));

	std::vector<double> values = {1, 2, 3, 2, 2, 3, 3, -1, -1, 4};
	CSRMatrix mat2(4, 4, positions, values);

	Vector v({1,2,3,0});
	Vector res = mat2 * v;
	ASSERT_EQ(mat2.numberOfRows(), res.getDimension());

	EXPECT_EQ(14, res[0]);
	EXPECT_EQ(6, res[1]);
	EXPECT_EQ(12, res[2]);
	EXPECT_EQ(-3, res[3]);

	// rectangular matrix
	positions = {std::make_pair(0,0), std::make_pair(1,2)};
	values = {1.0, 3.0};
	mat = CSRMatrix(2, 5, positions, values);

	v = {0,1,2,3,0};
	res = mat * v;

	ASSERT_EQ(2u, res.getDimension());
	EXPECT_EQ(0, res[0]);
	EXPECT_EQ(6, res[1]);
}

TEST(CSRMatrixGTest, testMatrixMultiplication) {
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values = {1, 2, 3, 2, 2, 3, 3, -1, -1, 4};

	positions.push_back(std::make_pair(0,0));
	positions.push_back(std::make_pair(0,1));
	positions.push_back(std::make_pair(0,2));
	positions.push_back(std::make_pair(1,0));
	positions.push_back(std::make_pair(1,1));
	positions.push_back(std::make_pair(2,0));
	positions.push_back(std::make_pair(2,2));
	positions.push_back(std::make_pair(2,3));
	positions.push_back(std::make_pair(3,2));
	positions.push_back(std::make_pair(3,3));

	//
	//				 1  2  3  0
	// 				 2  2  0  0
	// mat1 = mat2 = 3  0  3 -1
	//				 0  0 -1  4
	//
	CSRMatrix mat1(4, 4, positions, values);
	ASSERT_EQ(4u, mat1.numberOfRows());
	ASSERT_EQ(4u, mat1.numberOfColumns());

	CSRMatrix mat2(4, 4, positions, values);
	ASSERT_EQ(4u, mat2.numberOfRows());
	ASSERT_EQ(4u, mat2.numberOfColumns());

	//
	//			14  6  12  -3
	//			 6  8   6   0
	// result = 12  6  19  -7
	//			-3  0  -7  17
	//
	CSRMatrix result = mat1 * mat2;
	ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

	EXPECT_EQ(14, result(0,0));
	EXPECT_EQ(6, result(0,1));
	EXPECT_EQ(12, result(0,2));
	EXPECT_EQ(-3, result(0,3));
	EXPECT_EQ(6, result(1,0));
	EXPECT_EQ(8, result(1,1));
	EXPECT_EQ(6, result(1,2));
	EXPECT_EQ(0, result(1,3));
	EXPECT_EQ(12, result(2,0));
	EXPECT_EQ(6, result(2,1));
	EXPECT_EQ(19, result(2,2));
	EXPECT_EQ(-7, result(2,3));
	EXPECT_EQ(-3, result(3,0));
	EXPECT_EQ(0, result(3,1));
	EXPECT_EQ(-7, result(3,2));
	EXPECT_EQ(17, result(3,3));


	// rectangular matrices
	positions = {std::make_pair(0,0), std::make_pair(0,3), std::make_pair(1,2), std::make_pair(2,1), std::make_pair(2,3)};
	values = {1.0, 2.0, 1.0, 2.0, 4.0};
	mat1 = CSRMatrix(3, 4, positions, values);

	positions = {std::make_pair(0,0), std::make_pair(2,1), std::make_pair(3,0), std::make_pair(3,1)};
	values = {1.0, 0.5, 42.0, 1.0};
	mat2 = CSRMatrix(4, 2, positions, values);

	result = mat1 * mat2;

	EXPECT_EQ(85, result(0,0));
	EXPECT_EQ(2, result(0,1));
	EXPECT_EQ(0, result(1,0));
	EXPECT_EQ(0.5, result(1,1));
	EXPECT_EQ(168, result(2,0));
	EXPECT_EQ(4, result(2,1));
}

TEST(CSRMatrixGTest, testBigMatrixMultiplication) {
	METISGraphReader graphReader;
	Graph G = graphReader.read("input/PGPgiantcompo.graph");

	std::vector<std::pair<index,index>> positions;
	std::vector<double> values;

	G.forEdges([&](index i, index j, double value) {
		positions.push_back(std::make_pair(i,j));
		values.push_back(value);
	});

	CSRMatrix mat(G.upperNodeIdBound(), G.upperNodeIdBound(), positions, values);

	CSRMatrix result = mat * mat;
	ASSERT_EQ(mat.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat.numberOfColumns(), result.numberOfColumns());
}

TEST(CSRMatrixGTest, testTransposition) {
	//
	//	   1  2  3  1  1
	// 	   0  2  0  0  0
	// mat 4  0  3 -1  0
	//	   0  0  0  4 -1
	//
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values = {1, 2, 3, 1, 1, 2, 4, 3, -1, 4, -1};

	positions.push_back(std::make_pair(0,0));
	positions.push_back(std::make_pair(0,1));
	positions.push_back(std::make_pair(0,2));
	positions.push_back(std::make_pair(0,3));
	positions.push_back(std::make_pair(0,4));
	positions.push_back(std::make_pair(1,1));
	positions.push_back(std::make_pair(2,0));
	positions.push_back(std::make_pair(2,2));
	positions.push_back(std::make_pair(2,3));
	positions.push_back(std::make_pair(3,3));
	positions.push_back(std::make_pair(3,4));

	CSRMatrix mat(4, 5, positions, values);
	CSRMatrix matT = mat.transpose();

	EXPECT_EQ(5u, matT.numberOfRows());
	EXPECT_EQ(4u, matT.numberOfColumns());

	mat.forNonZeroElementsInRowOrder([&](index i, index j, double value) {
		EXPECT_EQ(value, matT(j,i));
	});
}

TEST(CSRMatrixGTest, testMatrixTransposeMatrixMultiplication) {
	std::vector<std::pair<index,index>> positions = {std::make_pair(0,0), std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,1), std::make_pair(2,0), std::make_pair(3,2)};
	std::vector<double> values = {1.0, 2.0, 3.0, 2.0, 3.0, -1.0};
	CSRMatrix A(4, 3, positions, values);

	positions = {std::make_pair(0,0), std::make_pair(1,0), std::make_pair(2,1), std::make_pair(3,1), std::make_pair(3,2)};
	values = {1.0, 3.0, -2.0, 5.0, -8.0};
	CSRMatrix B(4, 3, positions, values);

	CSRMatrix C = CSRMatrix::mTmMultiply(A, B);

	EXPECT_EQ(1, C(0,0));
	EXPECT_EQ(-6, C(0,1));
	EXPECT_EQ(0, C(0,2));
	EXPECT_EQ(8, C(1,0));
	EXPECT_EQ(0, C(1,1));
	EXPECT_EQ(0, C(1,2));
	EXPECT_EQ(3, C(2,0));
	EXPECT_EQ(-5, C(2,1));
	EXPECT_EQ(8, C(2,2));
}

TEST(CSRMatrixGTest, testMatrixMatrixTransposeMultiplication) {
	std::vector<std::pair<index,index>> positions = {std::make_pair(0,0), std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,1), std::make_pair(2,0), std::make_pair(3,2)};
	std::vector<double> values = {1.0, 2.0, 3.0, 2.0, 3.0, -1.0};
	CSRMatrix A(4, 3, positions, values);

	positions = {std::make_pair(0,0), std::make_pair(1,0), std::make_pair(2,1), std::make_pair(3,1), std::make_pair(3,2)};
	values = {1.0, 3.0, -2.0, 5.0, -8.0};
	CSRMatrix B(4, 3, positions, values);

	CSRMatrix C = CSRMatrix::mmTMultiply(A, B);

	EXPECT_EQ(1, C(0,0));
	EXPECT_EQ(3, C(0,1));
	EXPECT_EQ(-4, C(0,2));
	EXPECT_EQ(-14, C(0,3));
	EXPECT_EQ(-4, C(1,2));
	EXPECT_EQ(10, C(1,3));
	EXPECT_EQ(3, C(2,0));
	EXPECT_EQ(9, C(2,1));
	EXPECT_EQ(8, C(3,3));
	EXPECT_EQ(0, C(1,0));
	EXPECT_EQ(0, C(1,1));
	EXPECT_EQ(0, C(2,3));
}

TEST(CSRMatrixGTest, testMatrixTransposeVectorMultiplication) {
	std::vector<std::pair<index,index>> positions = {std::make_pair(0,0), std::make_pair(2,1)};
	std::vector<double> values = {1.0, 3.0};
	CSRMatrix mat(5, 2, positions, values);

	Vector v = {0,1,2,3,0};
	Vector res = CSRMatrix::mTvMultiply(mat, v);

	ASSERT_EQ(2u, res.getDimension());
	EXPECT_EQ(0, res[0]);
	EXPECT_EQ(6, res[1]);
}

} /* namespace NetworKit */
