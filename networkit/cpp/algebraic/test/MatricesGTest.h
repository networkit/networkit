/*
 * MatricesGTest.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_TEST_MATRICESGTEST_H_
#define NETWORKIT_CPP_ALGEBRAIC_TEST_MATRICESGTEST_H_

#include "gtest/gtest.h"
#include "../AlgebraicGlobals.h"
#include "../Vector.h"
#include "../../io/METISGraphReader.h"
#include "../MatrixTools.h"

namespace NetworKit {

class MatricesGTest : public testing::Test {
public:
	MatricesGTest() = default;
	virtual ~MatricesGTest() = default;

	template<class MATRIX>
	void testDimension();

	template<class MATRIX>
	void testNNZInRow();

	template<class MATRIX>
	void testRowAndColumnAccess();

	template<class MATRIX>
	void testDiagonalVector();

	template<class MATRIX>
	void testTranspose();

	template<class MATRIX>
	void testMatrixAddition();

	template<class MATRIX>
	void testMatrixSubtraction();

	template<class MATRIX>
	void testScalarMultiplication();

	template<class MATRIX>
	void testMatrixDivisionOperator();

	template<class MATRIX>
	void testMatrixVectorProduct();

	template<class MATRIX>
	void testMatrixMultiplication();

	template<class MATRIX>
	void testBigMatrixMultiplication();

	template<class MATRIX>
	void testLaplacianOfGraph();

	// TODO: Test other matrix classes

	// TODO: Test mmT multiplication, etc.!
};

template<class MATRIX>
void MatricesGTest::testDimension() {
	MATRIX mat(10);

	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = MATRIX(5, 10, 0.0);
	ASSERT_EQ(5u, mat.numberOfRows());
	ASSERT_EQ(10u, mat.numberOfColumns());

	mat = MATRIX(10, 5, 0.0);
	ASSERT_EQ(10u, mat.numberOfRows());
	ASSERT_EQ(5u, mat.numberOfColumns());
}

template<class MATRIX>
void MatricesGTest::testNNZInRow() {
	/*
	 * 1.0  0.0  2.0
	 * 4.0  0.0  0.0
	 * 0.0  0.0  0.0
	 * 0.0  0.0  2.0
	 */
	std::vector<Triplet> triplets = {{0,0,1.0}, {0,2,2.0}, {1,0,4.0}, {3,3,2.0}};

	MATRIX mat(4, triplets);
	EXPECT_EQ(2u, mat.nnzInRow(0));
	EXPECT_EQ(1u, mat.nnzInRow(1));
	EXPECT_EQ(0u, mat.nnzInRow(2));
	EXPECT_EQ(1u, mat.nnzInRow(3));
}

template<class MATRIX>
void MatricesGTest::testRowAndColumnAccess() {
	std::vector<Triplet> triplets;

	for (index i = 0; i < 1000; ++i) {
		triplets.push_back({3,i,(double)i});
	}

	triplets.push_back({10,10,42.123});

	MATRIX mat(1000, triplets);

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
	triplets.clear();
	triplets.push_back({4,9,11});
	triplets.push_back({0,0,42});

	mat = MATRIX(5, 10, triplets);
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

	triplets.clear();
	triplets.push_back({9,4,11});
	triplets.push_back({0,0,42});

	mat = MATRIX(10, 5, triplets);
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

template<class MATRIX>
void MatricesGTest::testDiagonalVector() {
	// 1  2  3  0
	// 2  2  0  0
	// 3  0  0 -1
	// 0  0 -1  4
	std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

	MATRIX mat(4,4,triplets);
	Vector diagonal = mat.diagonal();
	EXPECT_EQ(4u, diagonal.getDimension());
	EXPECT_EQ(1, diagonal[0]);
	EXPECT_EQ(2, diagonal[1]);
	EXPECT_EQ(0, diagonal[2]);
	EXPECT_EQ(4, diagonal[3]);

	// rectangular matrix
	//
	// 1 0 0 0 0
	// 0 0 3 0 0
	triplets = {{0,0,1}, {1,2,3}};
	mat = MATRIX(2,5,triplets);
	diagonal = mat.diagonal();
	EXPECT_EQ(2u, diagonal.getDimension());
	EXPECT_EQ(1, diagonal[0]);
	EXPECT_EQ(0, diagonal[1]);

	// rectangular matrix
	//
	// 1 0
	// 0 -3
	// 0 0
	// 0 0
	// 0 0
	triplets = {{0,0,1.0}, {1,1,-3.0}};
	mat = MATRIX(5,2,triplets);
	diagonal = mat.diagonal();
	EXPECT_EQ(2u, diagonal.getDimension());
	EXPECT_EQ(1, diagonal[0]);
	EXPECT_EQ(-3, diagonal[1]);
}

template<class MATRIX>
void MatricesGTest::testTranspose() {
	// 1  0  1  0
	// 2  2  0  0
	// 3  0  0 -1
	// 0  0  1  4
	std::vector<Triplet> triplets = {{0,0,1}, {0,2,1}, {1,0,2}, {1,1,2}, {2,0,3}, {2,3,-1}, {3,2,1}, {3,3,4}};

	MATRIX mat(4,4,triplets);
	MATRIX matT = mat.transpose();
	EXPECT_EQ(4u, matT.numberOfRows());
	EXPECT_EQ(4u, matT.numberOfColumns());

	EXPECT_EQ(1, matT(0,0));
	EXPECT_EQ(2, matT(0,1));
	EXPECT_EQ(3, matT(0,2));
	EXPECT_EQ(0, matT(0,3));
	EXPECT_EQ(0, matT(1,0));
	EXPECT_EQ(2, matT(1,1));
	EXPECT_EQ(0, matT(1,2));
	EXPECT_EQ(0, matT(1,3));
	EXPECT_EQ(1, matT(2,0));
	EXPECT_EQ(0, matT(2,1));
	EXPECT_EQ(0, matT(2,2));
	EXPECT_EQ(1, matT(2,3));
	EXPECT_EQ(0, matT(3,0));
	EXPECT_EQ(0, matT(3,1));
	EXPECT_EQ(-1, matT(3,2));
	EXPECT_EQ(4, matT(3,3));

	// rectangular matrix
	//
	// 1 0 0 0 0
	// 0 0 3 0 0
	triplets = {{0,0,1}, {1,2,3}};
	mat = MATRIX(2,5,triplets);
	matT = mat.transpose();
	EXPECT_EQ(5u, matT.numberOfRows());
	EXPECT_EQ(2u, matT.numberOfColumns());

	EXPECT_EQ(1, matT(0,0));
	EXPECT_EQ(0, matT(0,1));
	EXPECT_EQ(0, matT(1,0));
	EXPECT_EQ(0, matT(1,1));
	EXPECT_EQ(0, matT(2,0));
	EXPECT_EQ(3, matT(2,1));
	EXPECT_EQ(0, matT(3,0));
	EXPECT_EQ(0, matT(3,1));
	EXPECT_EQ(0, matT(4,0));
	EXPECT_EQ(0, matT(4,1));

	// rectangular matrix
	//
	// 1 0
	// 0 0
	// 0 3
	// 0 0
	// 0 0
	triplets = {{0,0,1.0}, {2,1,3.0}};
	mat = MATRIX(5,2,triplets);
	matT = mat.transpose();
	EXPECT_EQ(2u, matT.numberOfRows());
	EXPECT_EQ(5u, matT.numberOfColumns());

	EXPECT_EQ(1, matT(0,0));
	EXPECT_EQ(0, matT(0,1));
	EXPECT_EQ(0, matT(0,2));
	EXPECT_EQ(0, matT(0,3));
	EXPECT_EQ(0, matT(0,4));
	EXPECT_EQ(0, matT(1,0));
	EXPECT_EQ(0, matT(1,1));
	EXPECT_EQ(3, matT(1,2));
	EXPECT_EQ(0, matT(1,3));
	EXPECT_EQ(0, matT(1,4));
}

template<class MATRIX>
void MatricesGTest::testMatrixAddition() {
	std::vector<Triplet> triplets1;
	std::vector<Triplet> triplets2;

	for (index i = 0; i < 100; ++i) {
		triplets1.push_back({i,i,1});
		triplets2.push_back({i,i,(double)i});
	}

	triplets1.push_back({2,71, 1.8});
	triplets2.push_back({42,43,3.14});

	MATRIX mat1(100, triplets1);
	MATRIX mat2(100, triplets2);

	MATRIX result = mat1 + mat2;
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
	//
	// 1 0 0 0 0
	// 0 0 3 0 0
	std::vector<Triplet> triplets = {{0,0,1.0}, {1,2,3.0}};


	mat1 = MATRIX(2,5,triplets);

	// 0 0 1 0 0
	// 0 0 1 0 0
	triplets.clear();
	triplets = {{0,2,1.0}, {1,2,1.0}};

	mat2 = MATRIX(2,5,triplets);

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
	//
	// 1 0
	// 0 0
	// 0 3
	// 0 0
	// 0 0

	triplets.clear();
	triplets = {{0,0,1.0}, {2,1,3.0}};

	mat1 = MATRIX(5,2,triplets);

	// 0 0
	// 0 0
	// 1 1
	// 0 0
	// 0 0
	triplets.clear();
	triplets = {{2,0,1.0}, {2,1,1.0}};

	mat2 = MATRIX(5,2,triplets);

	result = mat1 + mat2;

	ASSERT_EQ(5u, result.numberOfRows());
	ASSERT_EQ(2u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(1, result(2,0));
	EXPECT_EQ(4, result(2,1));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(4,1));
}

template<class MATRIX>
void MatricesGTest::testMatrixSubtraction() {
	std::vector<Triplet> triplets1;
	std::vector<Triplet> triplets2;

	for (index i = 0; i < 100; ++i) {
		triplets1.push_back({i,i,1});
		triplets2.push_back({i,i,(double)i});
	}

	triplets1.push_back({2,71, 1.8});
	triplets2.push_back({42,43,3.14});

	MATRIX mat1(100, triplets1);
	MATRIX mat2(100, triplets2);

	MATRIX result = mat2 - mat1;
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
	//
	// 1 0 0 0 0
	// 0 0 3 0 0
	std::vector<Triplet> triplets = {{0,0,1.0}, {1,2,3.0}};


	mat1 = MATRIX(2,5,triplets);

	// 0 0 1 0 0
	// 0 0 1 0 0
	triplets.clear();
	triplets = {{0,2,1.0}, {1,2,1.0}};

	mat2 = MATRIX(2,5,triplets);

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
	//
	// 1 0
	// 0 0
	// 0 3
	// 0 0
	// 0 0

	triplets.clear();
	triplets = {{0,0,1.0}, {2,1,3.0}};

	mat1 = MATRIX(5,2,triplets);

	// 0 0
	// 0 0
	// 1 1
	// 0 0
	// 0 0
	triplets.clear();
	triplets = {{2,0,1.0}, {2,1,1.0}};

	mat2 = MATRIX(5,2,triplets);

	result = mat1 - mat2;

	ASSERT_EQ(5u, result.numberOfRows());
	ASSERT_EQ(2u, result.numberOfColumns());

	EXPECT_EQ(1, result(0,0));
	EXPECT_EQ(-1, result(2,0));
	EXPECT_EQ(2, result(2,1));

	EXPECT_EQ(0, result(0,1));
	EXPECT_EQ(0, result(4,1));
}

template<class MATRIX>
void MatricesGTest::testScalarMultiplication() {
	std::vector<Triplet> triplets;

	for (index i = 0; i < 100; ++i) {
		triplets.push_back({i,i,(double) i});
	}

	triplets.push_back({42,43,42.0});

	MATRIX mat(100, triplets);
	mat *= 2;
	ASSERT_EQ(100u, mat.numberOfRows());
	ASSERT_EQ(100u, mat.numberOfColumns());

	for (index i = 0; i < 100; ++i) {
		EXPECT_EQ(i*2, mat(i, i));
	}
	EXPECT_EQ(84.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 99));

	mat *= 0.5;

	for (index i = 0; i < 100; ++i) {
		EXPECT_EQ(i, mat(i, i));
	}
	EXPECT_EQ(42.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 99));

	// rectangular matrix
	//
	// 1 0 0 0 0
	// 0 0 3 0 0
	triplets = {{0,0,1.0}, {1,2,3.0}};
	mat = MATRIX(2,5,triplets);

	mat *= 2;

	EXPECT_EQ(2, mat(0,0));
	EXPECT_EQ(6, mat(1,2));
}

template<class MATRIX>
void MatricesGTest::testMatrixDivisionOperator() {
	std::vector<Triplet> triplets;

	for (index i = 0; i < 100; ++i) {
		triplets.push_back({i,i, (double) i});
	}

	triplets.push_back({42,43,42.0});

	MATRIX mat(100, triplets);
	mat /= (1.0 / 2.0);
	ASSERT_EQ(100u, mat.numberOfRows());
	ASSERT_EQ(100u, mat.numberOfColumns());

	for (index i = 0; i < 100; ++i) {
		EXPECT_EQ(i*2, mat(i, i));
	}
	EXPECT_EQ(84.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 99));

	mat /= 2;

	for (index i = 0; i < 100; ++i) {
		EXPECT_EQ(i, mat(i, i));
	}
	EXPECT_EQ(42.0, mat(42, 43));
	EXPECT_EQ(0.0, mat(55, 99));

	// rectangular matrix
	//
	// 1 0 0 0 0
	// 0 0 3 0 0
	triplets = {{0,0,1.0}, {1,2,3.0}};
	mat = MATRIX(2,5,triplets);

	mat /= 2;

	EXPECT_EQ(0.5, mat(0,0));
	EXPECT_EQ(1.5, mat(1,2));
}

template<class MATRIX>
void MatricesGTest::testMatrixVectorProduct() {
	std::vector<Triplet> triplets;

	for (index i = 0; i < 1000; ++i) {
		triplets.push_back({i,i, (double) i});
	}

	triplets.push_back({42,43,42.0});

	Vector vector(1000, 1.0);
	vector[500] = 3.5;

	MATRIX mat(1000, triplets);

	Vector result = mat * vector;
	ASSERT_EQ(mat.numberOfRows(), result.getDimension());

	for (index i = 0; i < 1000; ++i) {
		if (i != 500 && i != 42 && i != 43) {
			EXPECT_EQ(i, result[i]);
		}
	}

	EXPECT_EQ(42.0, mat(42, 43));
	EXPECT_EQ(84.0, result[42]);
	EXPECT_EQ(1750.0, result[500]);

	//		  1  2  3  0
	//        2  2  0  0
	// mat2 = 3  0  3 -1
	//		  0  0 -1  4
	triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

	MATRIX mat2(4, triplets);

	Vector v({1,2,3,0});
	Vector res = mat2 * v;
	ASSERT_EQ(mat2.numberOfRows(), res.getDimension());

	EXPECT_EQ(14, res[0]);
	EXPECT_EQ(6, res[1]);
	EXPECT_EQ(12, res[2]);
	EXPECT_EQ(-3, res[3]);

	// rectangular matrix
	//
	// 1 0 0 0 0
	// 0 0 3 0 0

	triplets = {{0,0,1}, {1,2,3}};
	mat = MATRIX(2,5,triplets);

	v = {0,1,2,3,0};
	res = mat * v;

	ASSERT_EQ(2u, res.getDimension());
	EXPECT_EQ(0, res[0]);
	EXPECT_EQ(6, res[1]);
}

template<class MATRIX>
void MatricesGTest::testMatrixMultiplication() {
	std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

	//
	//				 1  2  3  0
	// 				 2  2  0  0
	// mat1 = mat2 = 3  0  3 -1
	//				 0  0 -1  4
	//
	MATRIX mat1(4, triplets);
	ASSERT_EQ(4u, mat1.numberOfRows());
	ASSERT_EQ(4u, mat1.numberOfColumns());

	MATRIX mat2(4, triplets);
	ASSERT_EQ(4u, mat2.numberOfRows());
	ASSERT_EQ(4u, mat2.numberOfColumns());

	//
	//			14  6  12  -3
	//			 6  8   6   0
	// result = 12  6  19  -7
	//			-3  0  -7  17
	//
	MATRIX result = mat1 * mat2;
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
	//
	// 1 0 0 2
	// 0 0 1 0
	// 0 2 0 4
	triplets = {{0,0,1}, {0,3,2}, {1,2,1}, {2,1,2}, {2,3,4}};
	mat1 = MATRIX(3,4,triplets);

	//
	// 1  0
	// 0  0
	// 0  0.5
	// 42 1

	triplets = {{0,0,1}, {2,1,0.5}, {3,0,42}, {3,1,1}};
	mat2 = MATRIX(4,2, triplets);

	result = mat1 * mat2;
	ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat2.numberOfColumns(), result.numberOfColumns());

	EXPECT_EQ(85, result(0,0));
	EXPECT_EQ(2, result(0,1));
	EXPECT_EQ(0, result(1,0));
	EXPECT_EQ(0.5, result(1,1));
	EXPECT_EQ(168, result(2,0));
	EXPECT_EQ(4, result(2,1));
}

template<class MATRIX>
void MatricesGTest::testBigMatrixMultiplication() {
	METISGraphReader graphReader;
	MATRIX mat = MATRIX::adjacencyMatrix(graphReader.read("input/PGPgiantcompo.graph"));

	MATRIX result = mat * mat;
	ASSERT_EQ(mat.numberOfRows(), result.numberOfRows());
	ASSERT_EQ(mat.numberOfColumns(), result.numberOfColumns());
}


template<class MATRIX>
void MatricesGTest::testLaplacianOfGraph() {
	METISGraphReader graphReader;
	MATRIX mat = MATRIX::laplacianMatrix(graphReader.read("input/PGPgiantcompo.graph"));
	EXPECT_TRUE(MatrixTools::isLaplacian(mat));

	mat = MATRIX::laplacianMatrix(graphReader.read("input/power.graph"));
	EXPECT_TRUE(MatrixTools::isLaplacian(mat));
}


} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_TEST_MATRICESGTEST_H_ */
