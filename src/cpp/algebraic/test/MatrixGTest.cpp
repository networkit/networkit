/*
 * MatrixGTest.cpp
 *
 *  Created on: 16.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */
#ifndef NOGTEST

#include "MatrixGTest.h"

MatrixGTest::MatrixGTest() {
}

MatrixGTest::~MatrixGTest() {
}

TEST(MatrixGTest, tryMatrixDimension) {
	Matrix mat(10);

	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfRows());
}

TEST(MatrixGTest, tryRowAndColumnAccess) {
	std::vector<std::pair<int, int> > positions;
	std::vector<double> values;

	for (int i = 0; i < 1000; ++i) {
		positions.push_back(std::make_pair(3, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(10, 10));
	values.push_back(42.123);

	Matrix mat(1000, positions, values);

	Vector v = mat.row(3);
	for (int i = 0; i < 1000; ++i) {
		ASSERT_EQ(i, v(i));
	}

	v = mat.row(10);
	ASSERT_EQ(42.123, v(10));

	v = mat.column(10);
	ASSERT_EQ(10.0, v(3));
	ASSERT_EQ(42.123, v(10));

}

TEST(MatrixGTest, tryMatrixAddition) {
	std::vector<std::pair<NetworKit::node, NetworKit::node> > positions1;
	std::vector<std::pair<NetworKit::node, NetworKit::node> > positions2;
	std::vector<double> values1;
	std::vector<double> values2;

	for (int i = 0; i < 100000; ++i) {
		NetworKit::node n = static_cast<NetworKit::node>(i);
		positions1.push_back(std::make_pair(n, n));
		positions2.push_back(std::make_pair(n, n));
		values1.push_back(1);
		values2.push_back(i);
	}

	positions1.push_back(std::make_pair(static_cast<NetworKit::node>(2), static_cast<NetworKit::node>(71)));
	values1.push_back(1.8);

	positions2.push_back(std::make_pair(static_cast<NetworKit::node>(42), static_cast<NetworKit::node>(43)));
	values2.push_back(3.14);

	Matrix mat1(100000, positions1, values1);
	Matrix mat2(100000, positions2, values2);

	Matrix result = mat1 + mat2;
	ASSERT_EQ(0.0, result(10, 13));


	for (unsigned int i = 0; i < result.numberOfRows(); ++i) {
		ASSERT_EQ((i + 1), result(i, i));
	}
	ASSERT_EQ(1.8, result(2, 71));
	ASSERT_EQ(3.14, result(42, 43));

	ASSERT_EQ(0.0, result(3, 14));
}

TEST(MatrixGTest, tryMatrixSubtraction) {
	std::vector<std::pair<NetworKit::node, NetworKit::node> > positions1;
	std::vector<std::pair<NetworKit::node, NetworKit::node> > positions2;
	std::vector<double> values1;
	std::vector<double> values2;

	for (int i = 0; i < 100000; ++i) {
		NetworKit::node n = static_cast<NetworKit::node>(i);
		positions1.push_back(std::make_pair(n, n));
		positions2.push_back(std::make_pair(n, n));
		values1.push_back(1);
		values2.push_back(i);
	}

	positions1.push_back(std::make_pair(static_cast<NetworKit::node>(2), static_cast<NetworKit::node>(71)));
	values1.push_back(1.8);

	positions2.push_back(std::make_pair(static_cast<NetworKit::node>(42), static_cast<NetworKit::node>(43)));
	values2.push_back(3.14);

	Matrix mat1(100000, positions1, values1);
	Matrix mat2(100000, positions2, values2);

	Matrix result = mat2 - mat1;
	ASSERT_EQ(0.0, result(10, 13));


	for (unsigned int i = 0; i < result.numberOfRows(); ++i) {
		ASSERT_EQ(((int) i - 1), result(i, i));
	}
	ASSERT_EQ(-1.8, result(2, 71));
	ASSERT_EQ(3.14, result(42, 43));

	ASSERT_EQ(0.0, result(3, 14));
}

TEST(MatrixGTest, tryScalarMultiplication) {
	std::vector<std::pair<int, int> > positions;
	std::vector<double> values;

	for (int i = 0; i < 10000; ++i) {
		positions.push_back(std::make_pair(i, i));
		values.push_back(i);
	}

	positions.push_back(std::make_pair(42, 43));
	values.push_back(42.0);

	Matrix mat(10000, positions, values);
	mat *= 2;

	for (int i = 0; i < 10000; ++i) {
		ASSERT_EQ(i*2, mat(i, i));
	}
	ASSERT_EQ(84.0, mat(42, 43));
	ASSERT_EQ(0.0, mat(55, 199));

	mat *= 0.5;

	for (int i = 0; i < 10000; ++i) {
		ASSERT_EQ(i, mat(i, i));
	}
	ASSERT_EQ(42.0, mat(42, 43));
	ASSERT_EQ(0.0, mat(55, 199));
}

TEST(MatrixGTest, tryMatrixVectorProduct) {
	std::vector<std::pair<int, int> > mPositions;
	std::vector<double> mValues;

	for (int i = 0; i < 10000; ++i) {
		mPositions.push_back(std::make_pair(i, i));
		mValues.push_back(i);
	}

	mPositions.push_back(std::make_pair(42, 43));
	mValues.push_back(42.0);

	Vector vector(10000, 1.0);
	vector(500) = 3.5;

	Matrix mat(10000, mPositions, mValues);

	Vector result = mat * vector;

	for (int i = 0; i < 10000; ++i) {
		if (i != 500 && i != 42 && i != 43) {
			ASSERT_EQ(i, result(i));
		}
	}

	ASSERT_EQ(42.0, mat(42, 43));
	ASSERT_EQ(84.0, result(42));
	ASSERT_EQ(1750.0, result(500));


	std::vector<std::pair<int, int> > positions;
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
	ASSERT_EQ(14, res(0));
	ASSERT_EQ(6, res(1));
	ASSERT_EQ(12, res(2));
	ASSERT_EQ(-3, res(3));
}

TEST(MatrixGTest, tryMatrixMultiplication) {
	std::vector<std::pair<int, int> > positions;
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
	Matrix mat2(4, positions, values);

	//
	//			14  6  12  -3
	//			 6  8   6   0
	// result = 12  6  19  -7
	//			-3  0  -7  17
	//
	Matrix result = mat1 * mat2;

	ASSERT_EQ(14, result(0,0));
	ASSERT_EQ(6, result(0,1));
	ASSERT_EQ(12, result(0,2));
	ASSERT_EQ(-3, result(0,3));
	ASSERT_EQ(6, result(1,0));
	ASSERT_EQ(8, result(1,1));
	ASSERT_EQ(6, result(1,2));
	ASSERT_EQ(0, result(1,3));
	ASSERT_EQ(12, result(2,0));
	ASSERT_EQ(6, result(2,1));
	ASSERT_EQ(19, result(2,2));
	ASSERT_EQ(-7, result(2,3));
	ASSERT_EQ(-3, result(3,0));
	ASSERT_EQ(0, result(3,1));
	ASSERT_EQ(-7, result(3,2));
	ASSERT_EQ(17, result(3,3));
}



// special matrices tests
TEST(MatrixGTest, tryLaplacianMatrixConstruction) {
	NetworKit::Graph graph(6);
	graph.addEdge(0, 1);
	graph.addEdge(0, 4);
	graph.addEdge(1, 4);
	graph.addEdge(1, 2);
	graph.addEdge(2, 3);
	graph.addEdge(3, 4);
	graph.addEdge(3, 5);

	LaplacianMatrix laplaceMatrix(graph);
	ASSERT_EQ(2, laplaceMatrix(0,0));
	ASSERT_EQ(-1, laplaceMatrix(0,1));
	ASSERT_EQ(0, laplaceMatrix(0,2));
	ASSERT_EQ(0, laplaceMatrix(0,3));
	ASSERT_EQ(-1, laplaceMatrix(0,4));
	ASSERT_EQ(0, laplaceMatrix(0,5));
	ASSERT_EQ(3, laplaceMatrix(1,1));
	ASSERT_EQ(-1, laplaceMatrix(1,2));
	ASSERT_EQ(0, laplaceMatrix(1,3));
	ASSERT_EQ(-1, laplaceMatrix(1,4));
	ASSERT_EQ(0, laplaceMatrix(1,5));
	ASSERT_EQ(2, laplaceMatrix(2,2));
	ASSERT_EQ(-1, laplaceMatrix(2,3));
	ASSERT_EQ(0, laplaceMatrix(2,4));
	ASSERT_EQ(0, laplaceMatrix(2,5));
	ASSERT_EQ(3, laplaceMatrix(3,3));
	ASSERT_EQ(-1, laplaceMatrix(3,4));
	ASSERT_EQ(-1, laplaceMatrix(3,5));
	ASSERT_EQ(3, laplaceMatrix(4,4));
	ASSERT_EQ(0, laplaceMatrix(4,5));
	ASSERT_EQ(1, laplaceMatrix(5,5));
}

TEST(MatrixGTest, tryNormalizedLaplacianMatrixConstruction) {
	NetworKit::Graph graph(7);
	graph.addEdge(0, 1);
	graph.addEdge(0, 4);
	graph.addEdge(1, 4);
	graph.addEdge(1, 2);
	graph.addEdge(2, 3);
	graph.addEdge(3, 4);
	graph.addEdge(3, 5);

	NormalizedLaplacianMatrix normalizedLaplacianMatrix(graph);
	ASSERT_EQ(1, normalizedLaplacianMatrix(0,0));
	ASSERT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(0,1));
	ASSERT_EQ(0, normalizedLaplacianMatrix(0,2));
	ASSERT_EQ(0, normalizedLaplacianMatrix(0,3));
	ASSERT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(0,4));
	ASSERT_EQ(0, normalizedLaplacianMatrix(0,5));
	ASSERT_EQ(0, normalizedLaplacianMatrix(0,6));
	ASSERT_EQ(1, normalizedLaplacianMatrix(1,1));
	ASSERT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(1,2));
	ASSERT_EQ(0, normalizedLaplacianMatrix(1,3));
	ASSERT_EQ(-1.0 / 3.0, normalizedLaplacianMatrix(1,4));
	ASSERT_EQ(0, normalizedLaplacianMatrix(1,5));
	ASSERT_EQ(0, normalizedLaplacianMatrix(1,6));
	ASSERT_EQ(1, normalizedLaplacianMatrix(2,2));
	ASSERT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(2,3));
	ASSERT_EQ(0, normalizedLaplacianMatrix(2,4));
	ASSERT_EQ(0, normalizedLaplacianMatrix(2,5));
	ASSERT_EQ(0, normalizedLaplacianMatrix(2,6));
	ASSERT_EQ(1, normalizedLaplacianMatrix(3,3));
	ASSERT_EQ(-1.0 / 3.0, normalizedLaplacianMatrix(3,4));
	ASSERT_EQ(-1.0 / sqrt(3.0), normalizedLaplacianMatrix(3,5));
	ASSERT_EQ(0, normalizedLaplacianMatrix(3,6));
	ASSERT_EQ(1, normalizedLaplacianMatrix(4,4));
	ASSERT_EQ(0, normalizedLaplacianMatrix(4,5));
	ASSERT_EQ(0, normalizedLaplacianMatrix(4,6));
	ASSERT_EQ(1, normalizedLaplacianMatrix(5,5));
	ASSERT_EQ(0, normalizedLaplacianMatrix(5,6));
	ASSERT_EQ(0, normalizedLaplacianMatrix(6,6));
}

#endif

