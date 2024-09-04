/*
 * MatricesGTest.cpp
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <type_traits>

#include <gtest/gtest.h>

#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

namespace {

template <class Type>
class MatricesGTest : public testing::Test {
protected:
    using Matrix = Type;

    void SetUp() { graph = METISGraphReader{}.read("input/PGPgiantcompo.graph"); }

    static Matrix get4x4Matrix() {
        std::vector<Triplet> triplets = {{0, 0, 1}, {0, 1, 2},  {0, 2, 3},  {1, 0, 2}, {1, 1, 2},
                                         {2, 0, 3}, {2, 3, -1}, {3, 2, -1}, {3, 3, 4}};

        return {4, triplets};
    }

    Graph graph;
};

using Matrices = testing::Types<DynamicMatrix, CSRMatrix, DenseMatrix>;
TYPED_TEST_SUITE(MatricesGTest, Matrices, /*Comma needed for variadic macro.*/);

TYPED_TEST(MatricesGTest, testDimension) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat(10);

    EXPECT_EQ(mat.numberOfRows(), 10);
    EXPECT_EQ(mat.numberOfColumns(), 10);

    mat = Matrix(5, 10, 0.0);
    EXPECT_EQ(mat.numberOfRows(), 5);
    EXPECT_EQ(mat.numberOfColumns(), 10);

    mat = Matrix(10, 5, 0.0);
    EXPECT_EQ(mat.numberOfRows(), 10);
    EXPECT_EQ(mat.numberOfColumns(), 5);
}

TYPED_TEST(MatricesGTest, testNNZInRow) {
    /*
     * 1.0  0.0  2.0
     * 4.0  0.0  0.0
     * 0.0  0.0  0.0
     * 0.0  0.0  2.0
     */
    std::vector<Triplet> triplets = {{0, 0, 1.0}, {0, 2, 2.0}, {1, 0, 4.0}, {3, 3, 2.0}};

    typename TestFixture::Matrix mat(4, triplets);
    EXPECT_EQ(mat.nnzInRow(0), 2);
    EXPECT_EQ(mat.nnzInRow(1), 1);
    EXPECT_EQ(mat.nnzInRow(2), 0);
    EXPECT_EQ(mat.nnzInRow(3), 1);
}

TYPED_TEST(MatricesGTest, testRowAndColumnAccess) {
    const count num_triplets = 1000;
    std::vector<Triplet> triplets;
    triplets.reserve(num_triplets);

    for (index i = 0; i < num_triplets; ++i) {
        triplets.push_back({3, i, (double)i});
    }

    triplets.push_back({10, 10, 42.123});

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(num_triplets, triplets);

    Vector v = mat.row(3);
    ASSERT_EQ(mat.numberOfColumns(), v.getDimension());

    for (index i = 0; i < num_triplets; ++i) {
        EXPECT_EQ(i, v[i]);
    }

    v = mat.row(10);
    ASSERT_EQ(v.getDimension(), mat.numberOfColumns());
    ASSERT_TRUE(v.isTransposed());
    EXPECT_EQ(v[10], 42.123);

    v = mat.column(10);
    ASSERT_EQ(mat.numberOfRows(), v.getDimension());
    ASSERT_FALSE(v.isTransposed());

    EXPECT_EQ(v[3], 10.0);
    EXPECT_EQ(v[10], 42.123);

    // rectangular matrix
    // n x m (n < m)
    triplets.clear();
    triplets.push_back({4, 9, 11});
    triplets.push_back({0, 0, 42});

    mat = Matrix(5, 10, triplets);
    v = mat.row(0);
    ASSERT_EQ(v.getDimension(), 10);
    for (index i = 0; i < v.getDimension(); ++i) {
        if (i == 0) {
            EXPECT_EQ(v[i], 42);
        } else {
            EXPECT_EQ(v[i], 0);
        }
    }

    v = mat.column(9);
    ASSERT_EQ(v.getDimension(), 5);
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
    triplets.push_back({9, 4, 11});
    triplets.push_back({0, 0, 42});

    mat = Matrix(10, 5, triplets);
    v = mat.row(0);
    ASSERT_EQ(v.getDimension(), 5);
    for (index i = 0; i < v.getDimension(); ++i) {
        if (i == 0) {
            EXPECT_EQ(v[i], 42);
        } else {
            EXPECT_EQ(v[i], 0);
        }
    }

    v = mat.column(4);
    ASSERT_EQ(v.getDimension(), 10);
    for (index i = 0; i < v.getDimension(); ++i) {
        if (i == v.getDimension() - 1) {
            EXPECT_EQ(v[i], 11);
        } else {
            EXPECT_EQ(v[i], 0);
        }
    }
}

TYPED_TEST(MatricesGTest, testDiagonalVector) {
    // 1  2  3  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0 -1  4
    std::vector<Triplet> triplets = {{0, 0, 1}, {0, 1, 2},  {0, 2, 3},  {1, 0, 2}, {1, 1, 2},
                                     {2, 0, 3}, {2, 3, -1}, {3, 2, -1}, {3, 3, 4}};

    using Matrix = typename TestFixture::Matrix;

    Matrix mat(4, 4, triplets);
    Vector diagonal = mat.diagonal();
    EXPECT_EQ(diagonal.getDimension(), 4);
    EXPECT_EQ(diagonal[0], 1);
    EXPECT_EQ(diagonal[1], 2);
    EXPECT_EQ(diagonal[2], 0);
    EXPECT_EQ(diagonal[3], 4);

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0, 0, 1}, {1, 2, 3}};
    mat = Matrix(2, 5, triplets);
    diagonal = mat.diagonal();
    EXPECT_EQ(diagonal.getDimension(), 2);
    EXPECT_EQ(diagonal[0], 1);
    EXPECT_EQ(diagonal[1], 0);

    // rectangular matrix
    //
    // 1 0
    // 0 -3
    // 0 0
    // 0 0
    // 0 0
    triplets = {{0, 0, 1.0}, {1, 1, -3.0}};
    mat = Matrix(5, 2, triplets);
    diagonal = mat.diagonal();
    EXPECT_EQ(diagonal.getDimension(), 2);
    EXPECT_EQ(diagonal[0], 1);
    EXPECT_EQ(diagonal[1], -3);
}

TYPED_TEST(MatricesGTest, testTranspose) {
    // 1  0  1  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0  1  4
    std::vector<Triplet> triplets = {{0, 0, 1}, {0, 2, 1},  {1, 0, 2}, {1, 1, 2},
                                     {2, 0, 3}, {2, 3, -1}, {3, 2, 1}, {3, 3, 4}};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(4, 4, triplets);
    Matrix matT = mat.transpose();
    EXPECT_EQ(matT.numberOfRows(), 4);
    EXPECT_EQ(matT.numberOfColumns(), 4);

    EXPECT_EQ(matT(0, 0), 1);
    EXPECT_EQ(matT(0, 1), 2);
    EXPECT_EQ(matT(0, 2), 3);
    EXPECT_EQ(matT(0, 3), 0);
    EXPECT_EQ(matT(1, 0), 0);
    EXPECT_EQ(matT(1, 1), 2);
    EXPECT_EQ(matT(1, 2), 0);
    EXPECT_EQ(matT(1, 3), 0);
    EXPECT_EQ(matT(2, 0), 1);
    EXPECT_EQ(matT(2, 1), 0);
    EXPECT_EQ(matT(2, 2), 0);
    EXPECT_EQ(matT(2, 3), 1);
    EXPECT_EQ(matT(3, 0), 0);
    EXPECT_EQ(matT(3, 1), 0);
    EXPECT_EQ(matT(3, 2), -1);
    EXPECT_EQ(matT(3, 3), 4);

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0, 0, 1}, {1, 2, 3}};
    mat = Matrix(2, 5, triplets);
    matT = mat.transpose();
    EXPECT_EQ(matT.numberOfRows(), 5);
    EXPECT_EQ(matT.numberOfColumns(), 2);

    EXPECT_EQ(matT(0, 0), 1);
    EXPECT_EQ(matT(0, 1), 0);
    EXPECT_EQ(matT(1, 0), 0);
    EXPECT_EQ(matT(1, 1), 0);
    EXPECT_EQ(matT(2, 0), 0);
    EXPECT_EQ(matT(2, 1), 3);
    EXPECT_EQ(matT(3, 0), 0);
    EXPECT_EQ(matT(3, 1), 0);
    EXPECT_EQ(matT(4, 0), 0);
    EXPECT_EQ(matT(4, 1), 0);

    // rectangular matrix
    //
    // 1 0
    // 0 0
    // 0 3
    // 0 0
    // 0 0
    triplets = {{0, 0, 1.0}, {2, 1, 3.0}};
    mat = Matrix(5, 2, triplets);
    matT = mat.transpose();
    EXPECT_EQ(matT.numberOfRows(), 2);
    EXPECT_EQ(matT.numberOfColumns(), 5);

    EXPECT_EQ(matT(0, 0), 1);
    EXPECT_EQ(matT(0, 1), 0);
    EXPECT_EQ(matT(0, 2), 0);
    EXPECT_EQ(matT(0, 3), 0);
    EXPECT_EQ(matT(0, 4), 0);
    EXPECT_EQ(matT(1, 0), 0);
    EXPECT_EQ(matT(1, 1), 0);
    EXPECT_EQ(matT(1, 2), 3);
    EXPECT_EQ(matT(1, 3), 0);
    EXPECT_EQ(matT(1, 4), 0);
}

TYPED_TEST(MatricesGTest, testExtract) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = Matrix::adjacencyMatrix(this->graph);
    const index n = 500;
    std::vector<index> rows(n);
    std::vector<index> columns(n);

    for (index i = 0; i < n; ++i) {
        rows[i] = Aux::Random::integer(this->graph.numberOfNodes() - 1);
        columns[i] = Aux::Random::integer(this->graph.numberOfNodes() - 1);
    }

    Matrix subMatrix = mat.extract(rows, columns);
    ASSERT_EQ(rows.size(), subMatrix.numberOfRows());
    ASSERT_EQ(columns.size(), subMatrix.numberOfColumns());

    for (index i = 0; i < subMatrix.numberOfRows(); ++i) {
        for (index j = 0; j < subMatrix.numberOfColumns(); ++j) {
            EXPECT_EQ(mat(rows[i], columns[j]), subMatrix(i, j));
        }
    }

    // 1  0  1  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0  1  4
    std::vector<Triplet> triplets = {{0, 0, 1}, {0, 2, 1},  {1, 0, 2}, {1, 1, 2},
                                     {2, 0, 3}, {2, 3, -1}, {3, 2, 1}, {3, 3, 4}};

    mat = Matrix(4, 4, triplets);
    rows = {0, 2, 0};
    columns = {1, 0, 2, 1};

    //
    // 0 1 1 0
    // 0 3 0 0
    // 0 1 1 0
    //
    subMatrix = mat.extract(rows, columns);
    ASSERT_EQ(subMatrix.numberOfRows(), 3);
    ASSERT_EQ(subMatrix.numberOfColumns(), 4);

    EXPECT_EQ(subMatrix(0, 0), 0);
    EXPECT_EQ(subMatrix(0, 1), 1);
    EXPECT_EQ(subMatrix(0, 2), 1);
    EXPECT_EQ(subMatrix(0, 3), 0);
    EXPECT_EQ(subMatrix(1, 0), 0);
    EXPECT_EQ(subMatrix(1, 1), 3);
    EXPECT_EQ(subMatrix(1, 2), 0);
    EXPECT_EQ(subMatrix(1, 3), 0);
    EXPECT_EQ(subMatrix(2, 0), 0);
    EXPECT_EQ(subMatrix(2, 1), 1);
    EXPECT_EQ(subMatrix(2, 2), 1);
    EXPECT_EQ(subMatrix(2, 3), 0);
}

TYPED_TEST(MatricesGTest, testAssign) {
    // 1  0  1  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0  1  4
    std::vector<Triplet> triplets = {{0, 0, 1}, {0, 2, 1},  {1, 0, 2}, {1, 1, 2},
                                     {2, 0, 3}, {2, 3, -1}, {3, 2, 1}, {3, 3, 4}};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(4, 4, triplets);

    // -1 1
    //  3 0
    std::vector<Triplet> subTriplets = {{0, 0, -1}, {0, 1, 1}, {1, 0, 3}};
    Matrix sourceMat(2, 2, subTriplets);

    std::vector<index> rows = {2, 3};
    std::vector<index> columns = {0, 2};
    mat.assign(rows, columns, sourceMat);

    EXPECT_EQ(mat(0, 0), 1);
    EXPECT_EQ(mat(0, 1), 0);
    EXPECT_EQ(mat(0, 2), 1);
    EXPECT_EQ(mat(0, 3), 0);
    EXPECT_EQ(mat(1, 0), 2);
    EXPECT_EQ(mat(1, 1), 2);
    EXPECT_EQ(mat(1, 2), 0);
    EXPECT_EQ(mat(1, 3), 0);
    EXPECT_EQ(mat(2, 0), -1);
    EXPECT_EQ(mat(2, 1), 0);
    EXPECT_EQ(mat(2, 2), 1);
    EXPECT_EQ(mat(2, 3), -1);
    EXPECT_EQ(mat(3, 0), 3);
    EXPECT_EQ(mat(3, 1), 0);
    EXPECT_EQ(mat(3, 2), 0);
    EXPECT_EQ(mat(3, 3), 4);

    rows = {2, 3};
    columns = {2, 3};
    mat.assign(rows, columns, sourceMat);

    EXPECT_EQ(mat(0, 0), 1);
    EXPECT_EQ(mat(0, 1), 0);
    EXPECT_EQ(mat(0, 2), 1);
    EXPECT_EQ(mat(0, 3), 0);
    EXPECT_EQ(mat(1, 0), 2);
    EXPECT_EQ(mat(1, 1), 2);
    EXPECT_EQ(mat(1, 2), 0);
    EXPECT_EQ(mat(1, 3), 0);
    EXPECT_EQ(mat(2, 0), -1);
    EXPECT_EQ(mat(2, 1), 0);
    EXPECT_EQ(mat(2, 2), -1);
    EXPECT_EQ(mat(2, 3), 1);
    EXPECT_EQ(mat(3, 0), 3);
    EXPECT_EQ(mat(3, 1), 0);
    EXPECT_EQ(mat(3, 2), 3);
    EXPECT_EQ(mat(3, 3), 0);
}

TYPED_TEST(MatricesGTest, testApply) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = Matrix::adjacencyMatrix(this->graph);

    mat.apply([&](double value) { return 2 * value; });
    this->graph.forEdges([&](index i, index j, double value) { EXPECT_EQ(2 * value, mat(i, j)); });
}

TYPED_TEST(MatricesGTest, testMatrixAddition) {
    const int num_triplets = 100;
    std::vector<Triplet> triplets1(num_triplets);
    std::vector<Triplet> triplets2(num_triplets);

    for (index i = 0; i < num_triplets; ++i) {
        triplets1[i] = {i, i, 1};
        triplets2[i] = {i, i, (double)i};
    }

    triplets1.push_back({2, 71, 1.8});
    triplets2.push_back({42, 43, 3.14});

    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(num_triplets, triplets1);
    Matrix mat2(num_triplets, triplets2);

    Matrix result = mat1 + mat2;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(result(10, 13), 0.0);

    for (index i = 0; i < result.numberOfRows(); ++i) {
        EXPECT_EQ((i + 1), result(i, i));
    }
    EXPECT_EQ(result(2, 71), 1.8);
    EXPECT_EQ(result(42, 43), 3.14);

    EXPECT_EQ(result(3, 14), 0.0);

    // rectangular matrix
    // n x m (n < m)
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    std::vector<Triplet> triplets = {{0, 0, 1.0}, {1, 2, 3.0}};

    mat1 = Matrix(2, 5, triplets);

    // 0 0 1 0 0
    // 0 0 1 0 0
    triplets.clear();
    triplets = {{0, 2, 1.0}, {1, 2, 1.0}};

    mat2 = Matrix(2, 5, triplets);

    result = mat1 + mat2;

    ASSERT_EQ(result.numberOfRows(), 2);
    ASSERT_EQ(result.numberOfColumns(), 5);

    EXPECT_EQ(result(0, 0), 1);
    EXPECT_EQ(result(0, 2), 1);
    EXPECT_EQ(result(1, 2), 4);

    EXPECT_EQ(result(0, 1), 0);
    EXPECT_EQ(result(1, 4), 0);

    // rectangular matrix
    // n x m (n > m)
    //
    // 1 0
    // 0 0
    // 0 3
    // 0 0
    // 0 0

    triplets.clear();
    triplets = {{0, 0, 1.0}, {2, 1, 3.0}};

    mat1 = Matrix(5, 2, triplets);

    // 0 0
    // 0 0
    // 1 1
    // 0 0
    // 0 0
    triplets.clear();
    triplets = {{2, 0, 1.0}, {2, 1, 1.0}};

    mat2 = Matrix(5, 2, triplets);

    result = mat1 + mat2;

    ASSERT_EQ(result.numberOfRows(), 5);
    ASSERT_EQ(result.numberOfColumns(), 2);

    EXPECT_EQ(result(0, 0), 1);
    EXPECT_EQ(result(2, 0), 1);
    EXPECT_EQ(result(2, 1), 4);

    EXPECT_EQ(result(0, 1), 0);
    EXPECT_EQ(result(4, 1), 0);
}

TYPED_TEST(MatricesGTest, testMatrixAdditionEqual) {
    const int num_triplets = 100;
    std::vector<Triplet> triplets1(num_triplets);
    std::vector<Triplet> triplets2(num_triplets);

    for (index i = 0; i < num_triplets; ++i) {
        triplets1[i] = {i, i, 1};
        triplets2[i] = {i, i, static_cast<double>(i)};
    }

    triplets1.push_back({2, 71, 1.8});
    triplets2.push_back({42, 43, 3.14});

    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(num_triplets, triplets1);
    Matrix mat2(num_triplets, triplets2);

    Matrix result = mat1 + mat2;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(result(10, 13), 0.0);

    for (index i = 0; i < result.numberOfRows(); ++i) {
        EXPECT_EQ((i + 1), result(i, i));
    }
    EXPECT_EQ(result(2, 71), 1.8);
    EXPECT_EQ(result(42, 43), 3.14);

    EXPECT_EQ(result(3, 14), 0.0);

    // rectangular matrix
    // n x m (n < m)
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    std::vector<Triplet> triplets = {{0, 0, 1.0}, {1, 2, 3.0}};

    mat1 = Matrix(2, 5, triplets);

    // 0 0 1 0 0
    // 0 0 1 0 0
    triplets.clear();
    triplets = {{0, 2, 1.0}, {1, 2, 1.0}};

    mat2 = Matrix(2, 5, triplets);

    mat1 += mat2;

    ASSERT_EQ(mat1.numberOfRows(), 2);
    ASSERT_EQ(mat1.numberOfColumns(), 5);

    EXPECT_EQ(mat1(0, 0), 1);
    EXPECT_EQ(mat1(0, 2), 1);
    EXPECT_EQ(mat1(1, 2), 4);

    EXPECT_EQ(mat1(0, 1), 0);
    EXPECT_EQ(mat1(1, 4), 0);
}

TYPED_TEST(MatricesGTest, testMatrixSubtraction) {
    const int num_triplets = 100;
    std::vector<Triplet> triplets1(num_triplets);
    std::vector<Triplet> triplets2(num_triplets);

    for (index i = 0; i < num_triplets; ++i) {
        triplets1[i] = {i, i, 1};
        triplets2[i] = {i, i, (double)i};
    }

    triplets1.push_back({2, 71, 1.8});
    triplets2.push_back({42, 43, 3.14});

    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(num_triplets, triplets1);
    Matrix mat2(num_triplets, triplets2);

    Matrix result = mat2 - mat1;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(result(10, 13), 0.0);

    for (index i = 0; i < result.numberOfRows(); ++i) {
        EXPECT_EQ(((int)i - 1), result(i, i));
    }
    EXPECT_EQ(result(2, 71), -1.8);
    EXPECT_EQ(result(42, 43), 3.14);

    EXPECT_EQ(result(3, 14), 0.0);

    // rectangular matrix
    // n x m (n < m)
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    std::vector<Triplet> triplets = {{0, 0, 1.0}, {1, 2, 3.0}};

    mat1 = Matrix(2, 5, triplets);

    // 0 0 1 0 0
    // 0 0 1 0 0
    triplets.clear();
    triplets = {{0, 2, 1.0}, {1, 2, 1.0}};

    mat2 = Matrix(2, 5, triplets);

    result = mat1 - mat2;

    ASSERT_EQ(result.numberOfRows(), 2);
    ASSERT_EQ(result.numberOfColumns(), 5);

    EXPECT_EQ(result(0, 0), 1);
    EXPECT_EQ(result(0, 2), -1);
    EXPECT_EQ(result(1, 2), 2);

    EXPECT_EQ(result(0, 1), 0);
    EXPECT_EQ(result(1, 4), 0);

    // rectangular matrix
    // n x m (n > m)
    //
    // 1 0
    // 0 0
    // 0 3
    // 0 0
    // 0 0

    triplets.clear();
    triplets = {{0, 0, 1.0}, {2, 1, 3.0}};

    mat1 = Matrix(5, 2, triplets);

    // 0 0
    // 0 0
    // 1 1
    // 0 0
    // 0 0
    triplets.clear();
    triplets = {{2, 0, 1.0}, {2, 1, 1.0}};

    mat2 = Matrix(5, 2, triplets);

    result = mat1 - mat2;

    ASSERT_EQ(5u, result.numberOfRows());
    ASSERT_EQ(2u, result.numberOfColumns());

    EXPECT_EQ(result(0, 0), 1);
    EXPECT_EQ(result(2, 0), -1);
    EXPECT_EQ(result(2, 1), 2);

    EXPECT_EQ(result(0, 1), 0);
    EXPECT_EQ(result(4, 1), 0);
}

TYPED_TEST(MatricesGTest, testMatrixSubtractionEqual) {
    std::vector<Triplet> triplets1 = {{0, 0, 1.0}, {2, 1, 3.0}};
    std::vector<Triplet> triplets2 = {{2, 0, 1.0}, {2, 1, 1.0}};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(5, 2, triplets1);
    Matrix mat2(5, 2, triplets2);

    mat1 -= mat2;

    ASSERT_EQ(5u, mat1.numberOfRows());
    ASSERT_EQ(2u, mat1.numberOfColumns());

    EXPECT_EQ(mat1(0, 0), 1);
    EXPECT_EQ(mat1(2, 0), -1);
    EXPECT_EQ(mat1(2, 1), 2);

    EXPECT_EQ(mat1(0, 1), 0);
    EXPECT_EQ(mat1(4, 1), 0);
}

TYPED_TEST(MatricesGTest, testMatrixEquality) {
    std::vector<Triplet> triplets1 = {{0, 0, 1.0}};
    std::vector<Triplet> triplets2 = {{0, 0, 1.0}, {2, 1, 1.0}};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(3, 3, triplets1);
    Matrix mat2(3, 3, triplets2);

    EXPECT_FALSE(mat1 == mat2);

    mat1(2, 1) = 1.0; // add the missing value

    EXPECT_TRUE(mat1 == mat2);
}

TYPED_TEST(MatricesGTest, testMatrixAlmostEquality) {
    std::vector<Triplet> triplets1 = {{0, 0, 1.0}, {2, 1, 1.0}};
    std::vector<Triplet> triplets2 = {{0, 0, 1.0}, {2, 1, 1.1}};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(3, 3, triplets1);
    Matrix mat2(3, 3, triplets2);

    double eps = 0.01;
    EXPECT_FALSE(mat1.isApprox(mat2, eps)); // epsilon 0.01 < 0.1, therefore expecting false

    mat2(2, 1) = 1.001;

    EXPECT_TRUE(mat1.isApprox(mat2, eps)); // now epsilon 0.01 > 0.001
}

TYPED_TEST(MatricesGTest, testMatrixConversion) {
    std::vector<Triplet> triplets = {{2, 1, 1.0}, {3, 1, 2.0}, {0, 2, 1.5}, {0, 1, 1.0}};
    // ground truth
    DenseMatrix den(4, 4, triplets);
    CSRMatrix csr(4, 4, triplets);
    DynamicMatrix dyn(4, 4, triplets);

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(4, 4, triplets);

    DenseMatrix mat_to_den = MatrixTools::matrixTo<DenseMatrix>(mat);
    EXPECT_TRUE(den == mat_to_den);
    CSRMatrix mat_to_csr = MatrixTools::matrixTo<CSRMatrix>(mat);
    EXPECT_TRUE(csr == mat_to_csr);
    DynamicMatrix mat_to_dyn = MatrixTools::matrixTo<DynamicMatrix>(mat);
    EXPECT_TRUE(dyn == mat_to_dyn);
}

TYPED_TEST(MatricesGTest, testMatrixAssignmentOperator) {
    std::vector<Triplet> triplets = {{0, 0, 1.0}, {2, 1, 1.0}, {3, 1, 2.0}, {0, 2, 1.5}};
    std::vector<Triplet> triplets_plus_one = {
        {0, 0, 1.0}, {2, 1, 1.0}, {3, 1, 2.0}, {0, 2, 1.5}, {2, 3, 4.0}};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(4, 4, triplets);
    Matrix mat_plus_one(4, 4, triplets_plus_one);

    mat(2, 3) = 4;
    EXPECT_TRUE(mat == mat_plus_one);
}

TYPED_TEST(MatricesGTest, LaplacianMatrixRemoveNode) {
    using Matrix = typename TestFixture::Matrix;

    Graph G(4);
    G.addEdge(1, 2);
    auto L4 = Matrix::laplacianMatrix(G);
    Matrix expected4(4);
    expected4.setValue(1, 2, -1);
    expected4.setValue(2, 1, -1);
    expected4.setValue(1, 1, 1);
    expected4.setValue(2, 2, 1);

    EXPECT_TRUE(L4 == expected4);

    G.removeNode(3);
    auto L3 = Matrix::laplacianMatrix(G);

    Matrix expected3(3);
    expected3.setValue(1, 2, -1);
    expected3.setValue(2, 1, -1);
    expected3.setValue(1, 1, 1);
    expected3.setValue(2, 2, 1);
    EXPECT_TRUE(L3 == expected3);
}

TYPED_TEST(MatricesGTest, LaplacianMatrixRemoveEdge) {
    using Matrix = typename TestFixture::Matrix;

    Graph G(4);
    G.addEdge(1, 2);
    G.addEdge(0, 1);
    auto L = Matrix::laplacianMatrix(G);
    Matrix expected(4);
    expected.setValue(1, 2, -1);
    expected.setValue(2, 1, -1);
    expected.setValue(0, 1, -1);
    expected.setValue(1, 0, -1);
    expected.setValue(1, 1, 2);
    expected.setValue(0, 0, 1);
    expected.setValue(2, 2, 1);

    EXPECT_TRUE(L == expected);

    G.removeEdge(1, 2);
    auto LRemove = Matrix::laplacianMatrix(G);
    Matrix expectedRemove(4);
    expectedRemove.setValue(1, 0, -1);
    expectedRemove.setValue(0, 1, -1);
    expectedRemove.setValue(1, 1, 1);
    expectedRemove.setValue(0, 0, 1);
    EXPECT_TRUE(LRemove == expectedRemove);
}

TYPED_TEST(MatricesGTest, testScalarMultiplication) {
    const int num_triplets = 100;
    std::vector<Triplet> triplets(num_triplets);

    for (index i = 0; i < num_triplets; ++i) {
        triplets[i] = {i, i, (double)i};
    }

    triplets.push_back({42, 43, 42.0});

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(num_triplets, triplets);
    mat *= 2;
    ASSERT_EQ(mat.numberOfRows(), 100);
    ASSERT_EQ(mat.numberOfColumns(), 100);

    for (index i = 0; i < num_triplets; ++i) {
        EXPECT_EQ(i * 2, mat(i, i));
    }
    EXPECT_EQ(mat(42, 43), 84.0);
    EXPECT_EQ(mat(55, 99), 0.0);

    mat *= 0.5;

    for (index i = 0; i < num_triplets; ++i) {
        EXPECT_EQ(i, mat(i, i));
    }
    EXPECT_EQ(mat(42, 43), 42.0);
    EXPECT_EQ(mat(55, 99), 0.0);

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0, 0, 1.0}, {1, 2, 3.0}};
    mat = Matrix(2, 5, triplets);

    mat *= 2;

    EXPECT_EQ(mat(0, 0), 2);
    EXPECT_EQ(mat(1, 2), 6);
}

TYPED_TEST(MatricesGTest, testMatrixScalar) {
    std::vector<Triplet> triplets = {{0, 0, 1.0}, {2, 1, 3.0}};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(5, 2, triplets);
    // 1 0
    // 0 0
    // 0 3
    // 0 0
    // 0 0

    Matrix res = mat1 * 2;

    ASSERT_EQ(5u, res.numberOfRows());
    ASSERT_EQ(2u, res.numberOfColumns());

    EXPECT_EQ(res(0, 0), 2);
    EXPECT_EQ(res(2, 1), 6);
    EXPECT_EQ(res(0, 1), 0);
    EXPECT_EQ(res(4, 1), 0);
}

TYPED_TEST(MatricesGTest, testMatrixDivisionOperator) {
    const int num_triplets = 100;
    std::vector<Triplet> triplets(num_triplets);

    for (index i = 0; i < num_triplets; ++i) {
        triplets[i] = {i, i, (double)i};
    }

    triplets.push_back({42, 43, 42.0});

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(num_triplets, triplets);
    mat /= (1.0 / 2.0);
    ASSERT_EQ(mat.numberOfRows(), 100);
    ASSERT_EQ(mat.numberOfColumns(), 100);

    for (index i = 0; i < num_triplets; ++i) {
        EXPECT_EQ(i * 2, mat(i, i));
    }
    EXPECT_EQ(mat(42, 43), 84.0);
    EXPECT_EQ(mat(55, 99), 0.0);

    mat /= 2;

    for (index i = 0; i < num_triplets; ++i) {
        EXPECT_EQ(i, mat(i, i));
    }
    EXPECT_EQ(mat(42, 43), 42.0);
    EXPECT_EQ(mat(55, 99), 0.0);

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0, 0, 1.0}, {1, 2, 3.0}};
    mat = Matrix(2, 5, triplets);

    mat /= 2;

    EXPECT_EQ(mat(0, 0), 0.5);
    EXPECT_EQ(mat(1, 2), 1.5);
}

TYPED_TEST(MatricesGTest, testMatrixVectorProduct) {
    const int num_triplets = 1000;
    std::vector<Triplet> triplets(num_triplets);

    for (index i = 0; i < num_triplets; ++i) {
        triplets[i] = {i, i, (double)i};
    }

    triplets.push_back({42, 43, 42.0});

    Vector vector(num_triplets, 1.0);
    vector[500] = 3.5;

    using Matrix = typename TestFixture::Matrix;
    Matrix mat(num_triplets, triplets);

    Vector result = mat * vector;
    ASSERT_EQ(mat.numberOfRows(), result.getDimension());

    for (index i = 0; i < num_triplets; ++i) {
        if (i != 500 && i != 42 && i != 43) {
            EXPECT_EQ(i, result[i]);
        }
    }

    EXPECT_EQ(mat(42, 43), 42.0);
    EXPECT_EQ(result[42], 84.0);
    EXPECT_EQ(result[500], 1750.0);

    //		  1  2  3  0
    //        2  2  0  0
    // mat2 = 3  0  3 -1
    //		  0  0 -1  4
    triplets = {{0, 0, 1}, {0, 1, 2}, {0, 2, 3},  {1, 0, 2},  {1, 1, 2},
                {2, 0, 3}, {2, 2, 3}, {2, 3, -1}, {3, 2, -1}, {3, 3, 4}};

    Matrix mat2(4, triplets);

    Vector v({1, 2, 3, 0});
    Vector res = mat2 * v;
    ASSERT_EQ(mat2.numberOfRows(), res.getDimension());

    EXPECT_EQ(res[0], 14);
    EXPECT_EQ(res[1], 6);
    EXPECT_EQ(res[2], 12);
    EXPECT_EQ(res[3], -3);

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0

    triplets = {{0, 0, 1}, {1, 2, 3}};
    mat = Matrix(2, 5, triplets);

    v = {0, 1, 2, 3, 0};
    res = mat * v;

    ASSERT_EQ(res.getDimension(), 2);
    EXPECT_EQ(res[0], 0);
    EXPECT_EQ(res[1], 6);
}

TYPED_TEST(MatricesGTest, testMatrixMultiplication) {
    std::vector<Triplet> triplets = {{0, 0, 1}, {0, 1, 2}, {0, 2, 3},  {1, 0, 2},  {1, 1, 2},
                                     {2, 0, 3}, {2, 2, 3}, {2, 3, -1}, {3, 2, -1}, {3, 3, 4}};

    //
    //				 1  2  3  0
    // 				 2  2  0  0
    // mat1 = mat2 = 3  0  3 -1
    //				 0  0 -1  4
    //
    using Matrix = typename TestFixture::Matrix;
    Matrix mat1(4, triplets);
    ASSERT_EQ(mat1.numberOfRows(), 4);
    ASSERT_EQ(mat1.numberOfColumns(), 4);

    Matrix mat2(4, triplets);
    ASSERT_EQ(mat2.numberOfRows(), 4);
    ASSERT_EQ(mat2.numberOfColumns(), 4);

    //
    //			14  6  12  -3
    //			 6  8   6   0
    // result = 12  6  19  -7
    //			-3  0  -7  17
    //
    Matrix result = mat1 * mat2;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(result(0, 0), 14);
    EXPECT_EQ(result(0, 1), 6);
    EXPECT_EQ(result(0, 2), 12);
    EXPECT_EQ(result(0, 3), -3);
    EXPECT_EQ(result(1, 0), 6);
    EXPECT_EQ(result(1, 1), 8);
    EXPECT_EQ(result(1, 2), 6);
    EXPECT_EQ(result(1, 3), 0);
    EXPECT_EQ(result(2, 0), 12);
    EXPECT_EQ(result(2, 1), 6);
    EXPECT_EQ(result(2, 2), 19);
    EXPECT_EQ(result(2, 3), -7);
    EXPECT_EQ(result(3, 0), -3);
    EXPECT_EQ(result(3, 1), 0);
    EXPECT_EQ(result(3, 2), -7);
    EXPECT_EQ(result(3, 3), 17);

    // rectangular matrices
    //
    // 1 0 0 2
    // 0 0 1 0
    // 0 2 0 4
    triplets = {{0, 0, 1}, {0, 3, 2}, {1, 2, 1}, {2, 1, 2}, {2, 3, 4}};
    mat1 = Matrix(3, 4, triplets);

    //
    // 1  0
    // 0  0
    // 0  0.5
    // 42 1

    triplets = {{0, 0, 1}, {2, 1, 0.5}, {3, 0, 42}, {3, 1, 1}};
    mat2 = Matrix(4, 2, triplets);

    result = mat1 * mat2;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat2.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(result(0, 0), 85);
    EXPECT_EQ(result(0, 1), 2);
    EXPECT_EQ(result(1, 0), 0);
    EXPECT_EQ(result(1, 1), 0.5);
    EXPECT_EQ(result(2, 0), 168);
    EXPECT_EQ(result(2, 1), 4);
}

TYPED_TEST(MatricesGTest, testMatricesTransposedMultiplication) {
    using Matrix = typename TestFixture::Matrix;
    if (std::is_same<Matrix, DenseMatrix>::value) {
        GTEST_SKIP() << "Skipping mmT/mTm test for DenseMatrix";
    }
    std::vector<Triplet> triplets1 = {{0, 0, 1}, {0, 1, 2}, {0, 2, 3}, {1, 0, 2}, {1, 1, 2}};
    //           1  2  3
    // mat =     2  2  0
    //
    DynamicMatrix mat(2, 3, triplets1);
    ASSERT_EQ(mat.numberOfRows(), 2);
    ASSERT_EQ(mat.numberOfColumns(), 3);
    //           1  2
    // mat^T =   2  2
    //           3  0
    DynamicMatrix result = DynamicMatrix::mmTMultiply(mat, mat);
    // (mat*mat^T)
    //          14  6
    // result = 6   8
    EXPECT_EQ(result(0, 0), 14);
    EXPECT_EQ(result(0, 1), 6);
    EXPECT_EQ(result(1, 0), 6);
    EXPECT_EQ(result(1, 1), 8);

    result = DynamicMatrix::mTmMultiply(mat, mat);
    // (mat^T*mat)
    //          5  6  3
    // result = 6  8  6
    //          3  6  9
    EXPECT_EQ(result(0, 0), 5);
    EXPECT_EQ(result(0, 1), 6);
    EXPECT_EQ(result(0, 2), 3);
    EXPECT_EQ(result(1, 0), 6);
    EXPECT_EQ(result(1, 1), 8);
    EXPECT_EQ(result(1, 2), 6);
    EXPECT_EQ(result(2, 0), 3);
    EXPECT_EQ(result(2, 1), 6);
    EXPECT_EQ(result(2, 2), 9);
}

TYPED_TEST(MatricesGTest, testMatrixTransposedVectorMultiplication) {
    using Matrix = typename TestFixture::Matrix;
    if (std::is_same<Matrix, DenseMatrix>::value) {
        GTEST_SKIP() << "Skipping mmT/mTm test for DenseMatrix";
    }
    std::vector<Triplet> triplets1 = {{0, 0, 1}, {1, 0, 2}, {2, 0, 3}, {0, 1, 2}, {1, 1, 2}};
    DynamicMatrix mat(3, 2, triplets1);
    Vector vec({1, 2, 3});
    //
    //  mat^T    *    vec  = result
    //  1  2  3       1
    //  2  2  0  *    2    =  14
    //                3        6
    //
    Vector result = DynamicMatrix::mTvMultiply(mat, vec);
    EXPECT_EQ(result[0], 14);
    EXPECT_EQ(result[1], 6);
}

TYPED_TEST(MatricesGTest, testAdjacencyMatrix) {
    Graph G(6);
    G.addEdge(0, 0);
    G.addEdge(0, 1);
    G.addEdge(0, 4);
    G.addEdge(1, 2);
    G.addEdge(1, 4);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(3, 5);

    using Matrix = typename TestFixture::Matrix;
    Matrix mat = Matrix::adjacencyMatrix(G);

    // first row
    EXPECT_EQ(mat(0, 0), 1);
    EXPECT_EQ(mat(0, 1), 1);
    EXPECT_EQ(mat(0, 2), 0);
    EXPECT_EQ(mat(0, 3), 0);
    EXPECT_EQ(mat(0, 4), 1);
    EXPECT_EQ(mat(0, 5), 0);

    // third row
    EXPECT_EQ(mat(2, 0), 0);
    EXPECT_EQ(mat(2, 1), 1);
    EXPECT_EQ(mat(2, 2), 0);
    EXPECT_EQ(mat(2, 3), 1);
    EXPECT_EQ(mat(2, 4), 0);
    EXPECT_EQ(mat(2, 5), 0);

    // fifth row
    EXPECT_EQ(mat(4, 0), 1);
    EXPECT_EQ(mat(4, 1), 1);
    EXPECT_EQ(mat(4, 2), 0);
    EXPECT_EQ(mat(4, 3), 1);
    EXPECT_EQ(mat(4, 4), 0);
    EXPECT_EQ(mat(4, 5), 0);

    // directed, weighted G
    Graph dGraph(4, true, true);
    dGraph.addEdge(0, 1, 2);
    dGraph.addEdge(0, 0, 42);
    dGraph.addEdge(2, 3, -3);
    dGraph.addEdge(3, 2, 5);

    mat = Matrix::adjacencyMatrix(dGraph);
    ASSERT_EQ(dGraph.numberOfNodes(), mat.numberOfRows());
    ASSERT_EQ(dGraph.numberOfNodes(), mat.numberOfColumns());

    EXPECT_EQ(mat(0, 1), 2);
    EXPECT_EQ(mat(1, 0), 0);
    EXPECT_EQ(mat(0, 0), 42);
    EXPECT_EQ(mat(2, 3), -3);
    EXPECT_EQ(mat(3, 2), 5);

    // read lesmis G
    METISGraphReader graphReader;
    G = graphReader.read("input/lesmis.graph");

    // create AdjacencyMatrix
    mat = Matrix::adjacencyMatrix(G);

    G.forNodes([&](node u) {
        G.forNodes([&](node v) {
            if (G.hasEdge(u, v)) {
                EXPECT_EQ(G.weight(u, v), mat(u, v));
            } else {
                EXPECT_EQ(mat(u, v), 0.0);
            }
        });
    });
}

TYPED_TEST(MatricesGTest, testDiagonalMatrix) {
    Vector diagonal = {1, 0, 4, -1};

    using Matrix = typename TestFixture::Matrix;
    Matrix mat = Matrix::diagonalMatrix(diagonal);

    EXPECT_EQ(mat.numberOfRows(), 4);
    EXPECT_EQ(mat.numberOfColumns(), 4);

    EXPECT_EQ(mat(0, 0), 1);
    EXPECT_EQ(mat(1, 1), 0);
    EXPECT_EQ(mat(2, 2), 4);
    EXPECT_EQ(mat(3, 3), -1);

    for (index i = 0; i < mat.numberOfRows(); ++i) {
        for (index j = 0; j < mat.numberOfColumns(); ++j) {
            if (i != j) {
                EXPECT_EQ(mat(i, j), 0);
            }
        }
    }
}

TYPED_TEST(MatricesGTest, testIncidenceMatrix) {
    Graph G = Graph(5, true);
    G.addEdge(0, 1, 4.0);
    G.addEdge(0, 2, 9.0);
    G.addEdge(0, 3, 16.0);
    G.addEdge(2, 3, 1.0);
    G.addEdge(4, 1, 25.0);
    G.addEdge(4, 4, 1.0);

    G.indexEdges();

    using Matrix = typename TestFixture::Matrix;
    Matrix mat = Matrix::incidenceMatrix(G);
    ASSERT_EQ(G.numberOfNodes(), mat.numberOfRows());
    ASSERT_EQ(G.numberOfEdges(), mat.numberOfColumns());

    EXPECT_EQ(mat(0, 0), std::sqrt(G.weight(0, 1)));
    EXPECT_EQ(mat(1, 0), -std::sqrt(G.weight(0, 1)));
    for (uint64_t i = 2; i < mat.numberOfRows(); ++i) {
        EXPECT_EQ(mat(i, 0), 0.0);
    }

    EXPECT_EQ(mat(2, 1), -std::sqrt(G.weight(0, 2)));

    EXPECT_EQ(mat(3, 2), -std::sqrt(G.weight(0, 3)));
    EXPECT_EQ(mat(3, 3), -std::sqrt(G.weight(2, 3)));

    for (uint64_t i = 0; i < mat.numberOfRows(); ++i) {
        EXPECT_EQ(mat(i, 5), 0.0);
    }

    Vector row0 = mat.row(0);
    ASSERT_EQ(row0.getDimension(), mat.numberOfColumns());

    EXPECT_EQ(row0[0], std::sqrt(G.weight(0, 1)));
    EXPECT_EQ(row0[1], std::sqrt(G.weight(0, 2)));
    EXPECT_EQ(row0[2], std::sqrt(G.weight(0, 3)));
    for (uint64_t j = 3; j < row0.getDimension(); ++j) {
        EXPECT_EQ(row0[j], 0.0);
    }

    for (uint64_t j = 0; j < 5; ++j) {
        Vector column = mat.column(j);
        ASSERT_EQ(column.getDimension(), mat.numberOfRows());

        double sum = 0.0;
        for (uint64_t i = 0; i < column.getDimension(); ++i) {
            sum += column[i];
        }

        EXPECT_EQ(sum, 0.0);
    }

    Vector column5 = mat.column(5);
    ASSERT_EQ(column5.getDimension(), mat.numberOfRows());

    for (uint64_t i = 0; i < column5.getDimension(); ++i) {
        EXPECT_EQ(column5[i], 0.0);
    }

    Vector v = {12, 3, 9, 28, 0, -1};

    Vector result = mat * v;
    ASSERT_EQ(result.getDimension(), mat.numberOfRows());

    EXPECT_EQ(result[0], 69);
    EXPECT_EQ(result[1], -24);
    EXPECT_EQ(result[2], 19);
    EXPECT_EQ(result[3], -64);
    EXPECT_EQ(result[4], 0);
}

TYPED_TEST(MatricesGTest, testIncidenceMatrixDirected) {
    Graph G = Graph(5, true, true);
    G.addEdge(0, 1, 4.0);
    G.addEdge(0, 2, 9.0);
    G.addEdge(0, 3, 16.0);
    G.addEdge(2, 3, 1.0);
    G.addEdge(4, 1, 25.0);
    G.addEdge(4, 4, 1.0);

    G.indexEdges();

    using Matrix = typename TestFixture::Matrix;
    Matrix mat = Matrix::incidenceMatrix(G);
    ASSERT_EQ(G.numberOfNodes(), mat.numberOfRows());
    ASSERT_EQ(G.numberOfEdges(), mat.numberOfColumns());

    EXPECT_EQ(mat(0, 0), std::sqrt(G.weight(0, 1)));
    EXPECT_EQ(mat(1, 0), -std::sqrt(G.weight(0, 1)));
    for (uint64_t i = 2; i < mat.numberOfRows(); ++i) {
        EXPECT_EQ(mat(i, 0), 0.0);
    }

    EXPECT_EQ(mat(2, 1), -std::sqrt(G.weight(0, 2)));

    EXPECT_EQ(mat(3, 2), -std::sqrt(G.weight(0, 3)));
    EXPECT_EQ(mat(3, 3), -std::sqrt(G.weight(2, 3)));

    for (uint64_t i = 0; i < mat.numberOfRows(); ++i) {
        EXPECT_EQ(mat(i, 5), 0.0);
    }
}

TYPED_TEST(MatricesGTest, testLaplacianOfGraph) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = Matrix::laplacianMatrix(this->graph);

    EXPECT_TRUE(MatrixTools::isLaplacian(mat));

    mat = Matrix::laplacianMatrix(METISGraphReader{}.read("input/power.graph"));
    EXPECT_TRUE(MatrixTools::isLaplacian(mat));
}

TYPED_TEST(MatricesGTest, testNormalizedLaplacianOfGraph) {
    using Matrix = typename TestFixture::Matrix;
    if (std::is_same<Matrix, DenseMatrix>::value) {
        GTEST_SKIP() << "Skipping normalizedLaplacian test for DenseMatrix";
    }

    Graph G = METISGraphReader{}.read("input/power.graph");
    DynamicMatrix DynMat = DynamicMatrix::normalizedLaplacianMatrix(G);

    // check properties of normalizedLaplacian
    EXPECT_TRUE(MatrixTools::isSymmetric(DynMat));
    for (node u = 0; u < G.numberOfNodes(); ++u) {
        EXPECT_EQ(DynMat(u, u), 1);
    }
}

TYPED_TEST(MatricesGTest, testForElementsInRow) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    count nnz = 0;

    for (index row = 0; row < 4; ++row)
        mat.forElementsInRow(row, [&](index col, double value) {
            EXPECT_EQ(mat(row, col), value);
            if (value != mat.getZero())
                ++nnz;
        });

    EXPECT_EQ(nnz, mat.nnz());
}

TYPED_TEST(MatricesGTest, testForNonZeroElementsInRow) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    count nnz = 0;

    for (index row = 0; row < 4; ++row)
        mat.forNonZeroElementsInRow(row, [&](index col, double value) {
            EXPECT_EQ(mat(row, col), value);
            ++nnz;
        });

    EXPECT_EQ(nnz, mat.nnz());
}

TYPED_TEST(MatricesGTest, testParallelForElementsInRow) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    std::atomic<count> nnz{0};

    for (index row = 0; row < 4; ++row)
        mat.parallelForElementsInRow(row, [&](index col, double value) {
            EXPECT_EQ(mat(row, col), value);
            if (value != mat.getZero()) {
                nnz.fetch_add(1, std::memory_order_relaxed);
            }
        });

    EXPECT_EQ(mat.nnz(), nnz.load(std::memory_order_relaxed));
}

TYPED_TEST(MatricesGTest, testForElementsInRowOrder) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    count nnz = 0;

    mat.forElementsInRowOrder([&](index row, index col, double value) {
        EXPECT_EQ(mat(row, col), value);
        if (value != mat.getZero())
            ++nnz;
    });

    EXPECT_EQ(nnz, mat.nnz());
}

TYPED_TEST(MatricesGTest, testParallelForNonZeroElementsInRow) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    std::atomic<count> nnz{0};

    for (index row = 0; row < 4; ++row)
        mat.parallelForNonZeroElementsInRow(row, [&](index col, double value) {
            EXPECT_EQ(mat(row, col), value);
            nnz.fetch_add(1, std::memory_order_relaxed);
        });

    EXPECT_EQ(nnz.load(std::memory_order_relaxed), mat.nnz());
}

TYPED_TEST(MatricesGTest, testParallelForElementsInRowOrder) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    std::atomic<count> nnz{0};

    mat.parallelForElementsInRowOrder([&](index row, index col, double value) {
        EXPECT_EQ(mat(row, col), value);
        if (value != mat.getZero())
            nnz.fetch_add(1, std::memory_order_relaxed);
    });

    EXPECT_EQ(nnz.load(std::memory_order_relaxed), mat.nnz());
}

TYPED_TEST(MatricesGTest, testForNonZeroElementsInRowOrder) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    count nnz = 0;

    mat.forNonZeroElementsInRowOrder([&](index row, index col, double value) {
        EXPECT_EQ(mat(row, col), value);
        ++nnz;
    });

    EXPECT_EQ(nnz, mat.nnz());
}

TYPED_TEST(MatricesGTest, testParallelForNonZeroElementsInRowOrder) {
    using Matrix = typename TestFixture::Matrix;
    Matrix mat = this->get4x4Matrix();
    std::atomic<count> nnz{0};

    mat.parallelForNonZeroElementsInRowOrder([&](index row, index col, double value) {
        EXPECT_EQ(mat(row, col), value);
        nnz.fetch_add(1, std::memory_order_relaxed);
    });

    EXPECT_EQ(nnz.load(std::memory_order_relaxed), mat.nnz());
}

TYPED_TEST(MatricesGTest, testBigMatrixMultiplication) {
    using Matrix = typename TestFixture::Matrix;
    if (std::is_same<Matrix, DenseMatrix>::value) {
        GTEST_SKIP() << "Skipping big matrix multiplication test for DenseMatrix.";
    }

    Matrix mat = Matrix::adjacencyMatrix(this->graph);
    Matrix result = mat * mat;
    EXPECT_EQ(mat.numberOfRows(), result.numberOfRows());
    EXPECT_EQ(mat.numberOfColumns(), result.numberOfColumns());
}

TYPED_TEST(MatricesGTest, testMatrixToGraph) {
    using Matrix = typename TestFixture::Matrix;

    Matrix mat = Matrix::adjacencyMatrix(METISGraphReader{}.read("input/power.graph"));
    Graph G = MatrixTools::matrixToGraph(mat);

    EXPECT_EQ(G.numberOfNodes(), 4941);
    EXPECT_EQ(G.numberOfEdges(), 6594);
}

TYPED_TEST(MatricesGTest, testPrint) {
    Graph G(3);
    G.addEdge(1, 2);
    auto L = DynamicMatrix::laplacianMatrix(G);

    std::stringstream ss;

    ss << L;
    EXPECT_EQ(ss.str(), "0, 0, 0\n0, 1, -1\n0, -1, 1");
}

TYPED_TEST(MatricesGTest, test1by4) {
    using Matrix = typename TestFixture::Matrix;

    Matrix m1(4, 1, 0);
    Matrix m2(1, 4, 0);

    m1(0, 0) = 1;
    m1(1, 0) = 1;
    m1(2, 0) = 1;
    m1(3, 0) = 1;

    m2(0, 0) = 1;
    m2(0, 1) = 1;
    m2(0, 2) = 1;
    m2(0, 3) = 1;

    Matrix m3 = m1 * m2;
}

} // namespace

namespace {

using CSRMatrixGTest = testing::Test;

TEST_F(CSRMatrixGTest, testCSRMatrixSort) {
    /* 1 0 0 0
     * 2 3 0 0
     * 0 4 5 6
     * 7 8 9 10
     */
    std::vector<Triplet> triplets{{0, 0, 1},  {1, 1, 3}, {1, 0, 2}, // Insert unsorted columns
                                  {2, 2, 5},  {2, 1, 4}, {2, 3, 6}, {3, 1, 8},
                                  {3, 3, 10}, {3, 0, 7}, {3, 2, 9}};

    CSRGeneralMatrix<double> csr(4, 4, triplets);
    csr.sort();
}

} // namespace
} // namespace NetworKit
