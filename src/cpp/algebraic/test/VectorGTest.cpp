/*
 * VectorGTest.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */
#ifndef NOGTEST

#include "VectorGTest.h"

VectorGTest::VectorGTest() {
}

VectorGTest::~VectorGTest() {
}

TEST(VectorGTest, tryVectorConstruction) {
	Vector v1;
	ASSERT_EQ(0, v1.getDimension());

	Vector v2(10, 1.0);
	ASSERT_EQ(10, v2.getDimension());
	for (int i = 0; i < v2.getDimension(); i++) {
		ASSERT_EQ(1.0, v2(i));
	}

	std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector v3(values);
	ASSERT_EQ(5, v3.getDimension());
	for (int i = 0; i < v3.getDimension(); i++) {
		ASSERT_EQ(static_cast<double>((i + 1)), v3(i));
	}
}

TEST(VectorGTest, tryDimension) {
	std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector testVector(values);

	ASSERT_EQ(5, testVector.getDimension());
}

TEST(VectorGTest, tryAccessVectorElement) {
	std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector testVector(values);

	double third = testVector(2);
	ASSERT_EQ(3.0, third);

	testVector(2)++;
	third = testVector(2);
	ASSERT_EQ(4.0, third);
}

TEST(VectorGTest, tryInnerProduct) {
	std::vector<double> values1 = {1.0, 0.0, -1.0, -5.0, 2.0};
	Vector v1(values1);

	std::vector<double> values2 = {1.0, 2.0,  3.0,  0.0, 5.0};
	Vector v2(values2);

	double dotProduct = Vector::innerProduct(v1, v2);
	ASSERT_EQ(8.0, dotProduct);
}

TEST(VectorGTest, tryVectorComparisonOperators) {
	std::vector<double> values1 = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector v1(values1);

	std::vector<double> values2 = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector v2(values2);

	Vector v3(5, 1.0);

	ASSERT_TRUE(v1 == v2);
	ASSERT_FALSE(v1 != v2);

	ASSERT_FALSE(v1 == v3);
	ASSERT_TRUE(v1 != v3);
}

TEST(VectorGTest, tryVectorScalarMultiplication) {
	std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector testVector(values);

	Vector vectorScalar = testVector * 2.0;
	Vector scalarVector = 2.0 * testVector;

	for (int i = 0; i < testVector.getDimension(); i++) {
		ASSERT_EQ((i + 1) * 2.0, vectorScalar(i));
		ASSERT_EQ((i + 1) * 2.0, scalarVector(i));
	}
}

TEST(VectorGTest, tryVectorAddition) {
	std::vector<double> values1 = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector v1(values1);

	std::vector<double> values2 = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector v2(values2);

	Vector v3 = v1 + v2;
	ASSERT_EQ(5, v3.getDimension());
	for (int i = 0; i < v3.getDimension(); i++) {
		ASSERT_EQ(v1(i) + v2(i), v3(i));
	}

	v1 += v2;
	for (int i = 0; i < v1.getDimension(); i++) {
		ASSERT_EQ(v3(i), v1(i));
	}
}

TEST(VectorGTest, tryVectorSubtraction) {
	std::vector<double> values1 = {2.0, 4.0, 6.0, 8.0, 10.0};
	Vector v1(values1);

	std::vector<double> values2 = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector v2(values2);

	Vector v3 = v1 - v2;
	ASSERT_EQ(5, v3.getDimension());
	for (int i = 0; i < v3.getDimension(); i++) {
		ASSERT_EQ(v1(i) - v2(i), v3(i));
	}

	v1 -= v2;
	for (int i = 0; i < v1.getDimension(); i++) {
		ASSERT_EQ(v3(i), v1(i));
	}
}

TEST(VectorGTest, tryVectorIterators) {
	std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
	Vector v(values);

	// constant iterators
	int idx = 0;
	auto constantTester = [&](const double &element) {
		ASSERT_EQ(v(idx), element);
		idx++;
	};

	v.forElements(constantTester);

	auto constantParallelTester = [&](const int &idx, const double &element) {
		ASSERT_EQ(v(idx), element);
	};
	v.parallelForElements(constantParallelTester);

	// non-constant iterators
	auto nonConstantTester = [&](double &element) {
		element++;
	};

	v.forElements(nonConstantTester);
	for (int i = 0; i < v.getDimension(); i++) {
		ASSERT_EQ((i + 2.0), v(i));
	}

	auto nonConstantParallelTester = [](const int &idx, double &element) {
		element++;
	};

	v.parallelForElements(nonConstantParallelTester);
	for (int i = 0; i < v.getDimension(); i++) {
		ASSERT_EQ((i + 3.0), v(i));
	}
}


#endif
