/*
 * GaussSeidelRelaxationGTest.cpp
 *
 *  Created on: 03.11.2014
 *      Author: Michael
 */

#include "GaussSeidelRelaxationGTest.h"

namespace NetworKit {

TEST(GaussSeidelRelaxationGTest, trySolve) {
	std::vector<Vector> rows;
	rows.push_back({10, -1, 2, 0});
	rows.push_back({-1, 11, -1, 3});
	rows.push_back({2, -1, 10, -1});
	rows.push_back({0, 3, -1, 8});
	Matrix A(rows);

	Vector b = {6, 25, -11, 15};
	Vector x = {0, 0, 0, 0};

	GaussSeidelRelaxation solver;
	Vector result = solver.relax(A, b, x);

	EXPECT_EQ(1, std::round(result[0]));
	EXPECT_EQ(2, std::round(result[1]));
	EXPECT_EQ(-1, std::round(result[2]));
	EXPECT_EQ(1, std::round(result[3]));
}

TEST(GaussSeidelRelaxationGTest, tryIteration) {
	std::vector<Vector> rows;
	rows.push_back({10, -1, 2, 0});
	rows.push_back({-1, 11, -1, 3});
	rows.push_back({2, -1, 10, -1});
	rows.push_back({0, 3, -1, 8});
	Matrix A(rows);

	Vector b = {6, 25, -11, 15};
	Vector x = {0, 0, 0, 0};

	GaussSeidelRelaxation solver;
	Vector result = solver.relax(A, b, x, 1);

	EXPECT_TRUE(result[0] > 0);
	EXPECT_TRUE(result[1] > 1);
	EXPECT_TRUE(result[2] < 0);
	EXPECT_TRUE(result[3] > 0);
}

} /* namespace NetworKit */
