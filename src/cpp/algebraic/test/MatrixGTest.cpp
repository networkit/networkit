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

	ASSERT_EQ(10, mat.numberOfRows());
	ASSERT_EQ(10, mat.numberOfRows());
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
	ASSERT_EQ(0.0, result(10,13));


	for (int i = 0; i < result.numberOfRows(); ++i) {
		ASSERT_EQ((i + 1), result(i, i));
	}
	ASSERT_EQ(1.8, result(2, 71));
	ASSERT_EQ(3.14, result(42, 43));

	ASSERT_EQ(0.0, result (3, 14));
}

#endif

