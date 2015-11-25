/*
 * MatrixGTest.cpp
 *
 *  Created on: 16.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */
#ifndef NOGTEST

#include "MatrixGTest.h"

namespace NetworKit {

MatrixGTest::MatrixGTest() {
}

MatrixGTest::~MatrixGTest() {
}

TEST(MatrixGTest, testMatrixDimension) {
	Matrix mat(10);

	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = Matrix(5, 10);
	ASSERT_EQ(5u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = Matrix(10, 5);
	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(5u, mat.numberOfColumns());
}

TEST(MatrixGTest, testRowAndColumnAccess) {
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values;

	for (index i = 0; i < 1000; ++i) {
		positions.push_back(std::make_pair(3, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(10, 10));
	values.push_back(42.123);

	Matrix mat(1000, positions, values);

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

	mat = Matrix(5, 10, positions, values);
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

	mat = Matrix(10, 5, positions, values);
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

TEST(MatrixGTest, testMatrixAddition) {
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

	Matrix mat1(100000, positions1, values1);
	Matrix mat2(100000, positions2, values2);

	Matrix result = mat1 + mat2;
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
	std::vector<Vector> rows;
	rows.push_back({1,0,0,0,0});
	rows.push_back({0,0,3,0,0});

	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({0,0,1,0,0});
	rows.push_back({0,0,1,0,0});

	mat2 = Matrix(rows);

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
	rows.clear();
	rows.push_back({1,0});
	rows.push_back({0,0});
	rows.push_back({0,3});
	rows.push_back({0,0});
	rows.push_back({0,0});

	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({0,0});
	rows.push_back({0,0});
	rows.push_back({1,1});
	rows.push_back({0,0});
	rows.push_back({0,0});

	mat2 = Matrix(rows);

	result = mat1 + mat2;

	ASSERT_EQ(5u, result.numberOfRows());
	ASSERT_EQ(2u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(1, result(2,0));
	EXPECT_EQ(4, result(2,1));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(4,1));


	// non-matching dimensions
	//
	rows.clear();
	rows.push_back({1.0, 2.0});
	rows.push_back({2.0, 1.0});
	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({1.0, 0.0, 3.0});
	rows.push_back({0.0, 0.0, 1.0});
	rows.push_back({0.0, 0.0, 1.0});
	mat2 = Matrix(rows);

	EXPECT_THROW(mat1 + mat2, std::runtime_error);

	rows.clear();
	rows.push_back({1.0, 2.0, 0.0});
	rows.push_back({2.0, 1.0, 0.0});
	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({1.0, 0.0, 3.0});
	rows.push_back({0.0, 0.0, 1.0});
	rows.push_back({0.0, 0.0, 1.0});
	mat2 = Matrix(rows);

	EXPECT_THROW(mat1 + mat2, std::runtime_error);

	rows.clear();
	rows.push_back({1.0, 2.0});
	rows.push_back({2.0, 1.0});
	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({1.0, 0.0, 3.0});
	rows.push_back({0.0, 0.0, 1.0});
	mat2 = Matrix(rows);

	EXPECT_THROW(mat1 + mat2, std::runtime_error);
}

TEST(MatrixGTest, testMatrixSubtraction) {
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

	Matrix mat1(100000, positions1, values1);
	Matrix mat2(100000, positions2, values2);

	Matrix result = mat2 - mat1;
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
	std::vector<Vector> rows;
	rows.push_back({1,0,0,0,0});
	rows.push_back({0,0,3,0,0});

	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({0,0,1,0,0});
	rows.push_back({0,0,1,0,0});

	mat2 = Matrix(rows);

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
	rows.clear();
	rows.push_back({1,0});
	rows.push_back({0,0});
	rows.push_back({0,3});
	rows.push_back({0,0});
	rows.push_back({0,0});

	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({0,0});
	rows.push_back({0,0});
	rows.push_back({1,1});
	rows.push_back({0,0});
	rows.push_back({0,0});

	mat2 = Matrix(rows);

	result = mat1 - mat2;

	ASSERT_EQ(5u, result.numberOfRows());
	ASSERT_EQ(2u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(-1, result(2,0));
	EXPECT_EQ(2, result(2,1));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(4,1));


	// non-matching dimensions
	//
	rows.clear();
	rows.push_back({1.0, 2.0});
	rows.push_back({2.0, 1.0});
	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({1.0, 0.0, 3.0});
	rows.push_back({0.0, 0.0, 1.0});
	rows.push_back({0.0, 0.0, 1.0});
	mat2 = Matrix(rows);

	EXPECT_THROW(mat1 - mat2, std::runtime_error);

	rows.clear();
	rows.push_back({1.0, 2.0, 0.0});
	rows.push_back({2.0, 1.0, 0.0});
	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({1.0, 0.0, 3.0});
	rows.push_back({0.0, 0.0, 1.0});
	rows.push_back({0.0, 0.0, 1.0});
	mat2 = Matrix(rows);

	EXPECT_THROW(mat1 - mat2, std::runtime_error);

	rows.clear();
	rows.push_back({1.0, 2.0});
	rows.push_back({2.0, 1.0});
	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({1.0, 0.0, 3.0});
	rows.push_back({0.0, 0.0, 1.0});
	mat2 = Matrix(rows);

	EXPECT_THROW(mat1 - mat2, std::runtime_error);
}

TEST(MatrixGTest, testScalarMultiplication) {
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values;

	for (index i = 0; i < 10000; ++i) {
		positions.push_back(std::make_pair(i, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(42, 43));
	values.push_back(42.0);

	Matrix mat(10000, positions, values);
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
	std::vector<Vector> rows;
	rows.push_back({1,0,0,0,0});
	rows.push_back({0,0,3,0,0});
	mat = Matrix(rows);

	mat *= 2;

	EXPECT_EQ(2, mat(0,0));
	EXPECT_EQ(6, mat(1,2));
}

TEST(MatrixGTest, testMatrixDivisionOperator) {
	std::vector<std::pair<index, index> > positions;
	std::vector<double> values;

	for (index i = 0; i < 10000; ++i) {
		positions.push_back(std::make_pair(i, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(42, 43));
	values.push_back(42.0);

	Matrix mat(10000, positions, values);
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
	std::vector<Vector> rows;
	rows.push_back({1,0,0,0,0});
	rows.push_back({0,0,3,0,0});
	mat = Matrix(rows);

	mat /= 2;

	EXPECT_EQ(0.5, mat(0,0));
	EXPECT_EQ(1.5, mat(1,2));
}

TEST(MatrixGTest, testMatrixVectorProduct) {
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

	Matrix mat(10000, mPositions, mValues);

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
	Matrix mat2(4, positions, values);

	Vector v({1,2,3,0});
	Vector res = mat2 * v;
	ASSERT_EQ(mat2.numberOfRows(), res.getDimension());

	EXPECT_EQ(14, res[0]);
	EXPECT_EQ(6, res[1]);
	EXPECT_EQ(12, res[2]);
	EXPECT_EQ(-3, res[3]);

	// rectangular matrix
	std::vector<Vector> rows;
	rows.push_back({1,0,0,0,0});
	rows.push_back({0,0,3,0,0});
	mat = Matrix(rows);

	v = {0,1,2,3,0};
	res = mat * v;

	ASSERT_EQ(2u, res.getDimension());
	EXPECT_EQ(0, res[0]);
	EXPECT_EQ(6, res[1]);

	EXPECT_THROW(mat * v.transpose(), std::runtime_error);

	Vector v1 = {1.0, 2.0};
	EXPECT_THROW(mat * v1, std::runtime_error);
}

TEST(MatrixGTest, testMatrixMultiplication) {
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
	Matrix mat1(4, positions, values);
	ASSERT_EQ(4u, mat1.numberOfRows());
	ASSERT_EQ(4u, mat1.numberOfColumns());

	Matrix mat2(4, positions, values);
	ASSERT_EQ(4u, mat2.numberOfRows());
	ASSERT_EQ(4u, mat2.numberOfColumns());

	//
	//			14  6  12  -3
	//			 6  8   6   0
	// result = 12  6  19  -7
	//			-3  0  -7  17
	//
	Matrix result = mat1 * mat2;
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
	std::vector<Vector> rows;
	rows.push_back({1,0,0,2});
	rows.push_back({0,0,1,0});
	rows.push_back({0,2,0,4});

	mat1 = Matrix(rows);

	rows.clear();
	rows.push_back({1,0});
	rows.push_back({0,0});
	rows.push_back({0,0.5});
	rows.push_back({42,1});

	mat2 = Matrix(rows);

	result = mat1 * mat2;

	EXPECT_EQ(85, result(0,0));
	EXPECT_EQ(2, result(0,1));
	EXPECT_EQ(0, result(1,0));
	EXPECT_EQ(0.5, result(1,1));
	EXPECT_EQ(168, result(2,0));
	EXPECT_EQ(4, result(2,1));


	// non-matching dimensions
	//
	rows.clear();
	rows.push_back({1, 0});
	rows.push_back({0, 0});
	rows.push_back({0, 0.5});
	mat2 = Matrix(rows);

	EXPECT_THROW(mat1 * mat2, std::runtime_error);
}

TEST(MatrixGTest, testBigMatrixMultiplication) {
	METISGraphReader graphReader;
	AdjacencyMatrix mat(graphReader.read("input/PGPgiantcompo.graph"));

	Matrix result = mat * mat;
	ASSERT_EQ(mat.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat.numberOfColumns(), result.numberOfColumns());
}


} /* namespace NetworKit */

#endif

