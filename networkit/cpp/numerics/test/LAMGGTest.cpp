/*
 * LAMGGTest.cpp
 *
 *  Created on: 24.04.2024
 *      Author: Lukas
 */

#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>

namespace NetworKit {

// list of variants:
// - solver function: solve, parallelSolve
// - setup function: setupConnected, matrix, matrix graph, matrix graph decomp, graph decomp, graph
// - components of G: one component, multiple components
// - parallelism: single thread, multiple threads

class LAMGGTest
    // params: Solver, Setup, Components, Parallelism
    : public testing::TestWithParam<
          std::tuple<std::string, std::string, std::string, std::string>> {
protected:
    // returns a graph and a vector of RHSs
    std::tuple<Graph, std::vector<Vector>> testData() const;
    template <class Matrix>
    void runTest();
};

INSTANTIATE_TEST_SUITE_P(LAMGGTest, LAMGGTest,
                         testing::Combine(testing::Values("solve", "parallelSolve", "loaderSolve"),
                                          testing::Values("setupConnected", "matrix",
                                                          "matrix graph", "matrix graph decomp",
                                                          "graph decomp", "graph"),
                                          testing::Values("one", "two"),
                                          testing::Values("single", "four")));

inline bool vector_almost_equal(const Vector &lhs, const Vector &rhs) {
    if (lhs.getDimension() != rhs.getDimension())
        return false;
    for (size_t i = 0; i < lhs.getDimension(); ++i) {
        if (std::isnan(lhs[i]) || std::isnan(rhs[i]))
            return false;
        if (std::abs(lhs[i] - rhs[i]) > 1e-4)
            return false;
    }
    return true;
}

std::tuple<Graph, std::vector<Vector>> LAMGGTest::testData() const {
    auto [solveFn, setupFn, components, parallelism] = GetParam();

    Graph G(6);
    std::vector<Vector> rhss;
    std::vector<Vector> solutions;
    if (components == "one") {
        /*
        0 - 1 - 2
        |       |
        3 - 4 - 5

        L:
         2 -1  0 -1  0  0
        -1  2 -1  0  0  0
         0 -1  2  0  0 -1
        -1  0  0  2 -1  0
         0  0  0 -1  2 -1
         0  0 -1  0 -1  2
        */

        G.addEdge(0, 1);
        G.addEdge(0, 3);
        G.addEdge(1, 2);
        G.addEdge(2, 5);
        G.addEdge(3, 4);
        G.addEdge(4, 5);

        rhss = {
            {-4, 0, -2, 2, 0, 4},  {-8, 2, 2, 14, -8, -10}, {4, 0, -2, 2, 0, 4},
            {8, 2, 10, 0, -8, -1}, {4, 0, 6, 9, 0, 4},      {8, 8, -10, 5, -9, -3},
            {4, 0, 7, 2, 0, -4},   {8, 9, -7, 14, 8, -2},
        };

        // make sure that each vector sums up to 0
        for (auto &v : rhss) {
            double sum = 0;
            v.forElements([&](double value) { sum += value; });
            v[0] -= sum;
        }

        return {G, rhss};
    }

    if (components == "two") {
        /*
       0 - 1 - 2

       3 - 4 - 5

       L:
        1 -1  0  0  0  0
       -1  2 -1  0  0  0
        0 -1  1  0  0  0
        0  0  0  1 -1  0
        0  0  0 -1  2 -1
        0  0  0  0 -1  1
       */

        G.addEdge(0, 1);
        G.addEdge(1, 2);
        G.addEdge(3, 4);
        G.addEdge(4, 5);

        rhss = {
            {-1, 3, -2, -1, -3, 4},    {0, 0, 0, -1, -3, 4},       {-1, 3, -2, 0, 0, 0},
            {-10, 0, 10, 16, -6, -10}, {1, 3, 2, 1, 3, 4},         {0, 0, 0, -1, 3, 4},
            {1, 3, -2, 0, 0, 0},       {10, 0, -10, -16, -6, -10},
        };

        // make sure that each component sums up to 0
        for (auto &v : rhss) {
            v[0] -= v[0] + v[1] + v[2];
            v[3] -= v[3] + v[4] + v[5];
        }

        return {G, rhss};
    }
    throw std::logic_error("unknown number of components.");
}

TEST_P(LAMGGTest, testLamgCSRMatrixVariants) {
    runTest<CSRMatrix>();
}
TEST_P(LAMGGTest, testLamgDenseMatrixVariants) {
    runTest<DenseMatrix>();
}
TEST_P(LAMGGTest, testLamgDynamicMatrixVariants) {
    runTest<DynamicMatrix>();
}

template <class Matrix>
void LAMGGTest::runTest() {
    auto [solveFn, setupFn, components, parallelism] = GetParam();

    if (parallelism == "single")
        omp_set_num_threads(1);
    else if (parallelism == "four")
        omp_set_num_threads(4);
    else
        throw std::logic_error("unhandled variant!");

    // skip combinations that do not work
    if (setupFn == "setupConnected" && components != "one")
        return;

    // auto [G, rhss] = testData();
    auto data = testData();
    auto G = std::get<0>(data);
    auto rhss = std::get<1>(data);

    auto L = Matrix::laplacianMatrix(G);

    Lamg<Matrix> lamg;

    ParallelConnectedComponents pcc(G);
    pcc.run();

    if (setupFn == "setupConnected")
        lamg.setupConnected(L);
    else if (setupFn == "matrix")
        lamg.setup(L);
    else if (setupFn == "matrix graph")
        lamg.setup(L, G);
    else if (setupFn == "matrix graph decomp")
        lamg.setup(L, G, pcc);
    else if (setupFn == "graph decomp")
        lamg.setup(G, pcc);
    else if (setupFn == "graph")
        lamg.setup(G);
    else
        throw std::logic_error("unhandled variant!");

    // for loader solve variant
    std::atomic<int> numLoaderCalls = 0;
    const auto rhsLoader = [&rhss, &numLoaderCalls](count i, Vector &rhs) -> Vector & {
        rhs = rhss[i];
        ++numLoaderCalls;
        return rhs;
    };

    std::atomic<int> numProcessorCalls = 0;
    const auto resultProcessor = [&](count i, const Vector &result) {
        EXPECT_TRUE(vector_almost_equal(L * result, rhss[i]))
            << "Lamg result: " << result << "\nL * result: " << L * result << "\nrhs: " << rhss[i];
        ++numProcessorCalls;
    };

    if (solveFn == "solve") {
        for (size_t i = 0; i < rhss.size(); ++i) {
            const auto &rhs = rhss[i];

            Vector result(rhs.getDimension());

            lamg.solve(rhs, result);

            EXPECT_TRUE(vector_almost_equal(L * result, rhs))
                << "Lamg result: " << result << "\nL * result: " << L * result << "\nrhs: " << rhs;
        }
    } else if (solveFn == "parallelSolve") {
        std::vector<Vector> results(rhss.size(), Vector(6));
        lamg.parallelSolve(rhss, results);
        for (size_t i = 0; i < rhss.size(); ++i) {
            const auto &rhs = rhss[i];
            const auto &result = results[i];
            EXPECT_TRUE(vector_almost_equal(L * result, rhs))
                << "Lamg result: " << result << "\nL * result: " << L * result << "\nrhs: " << rhs;
        }
    } else if (solveFn == "loaderSolve") {
        lamg.parallelSolve(rhsLoader, resultProcessor, {rhss.size(), L.numberOfColumns()});
        EXPECT_EQ(numProcessorCalls, rhss.size());
        EXPECT_EQ(numLoaderCalls, rhss.size());
    } else
        throw std::logic_error("unhandled variant!");
}

TEST(LAMGGTest, starEliminatedToOneNode) {
    omp_set_num_threads(1);

    // A star eliminates every leaf and leaves a one-row Laplacian.
    Graph graph(201);
    for (node v = 1; v < 201; ++v)
        graph.addEdge(0, v);

    const auto L = CSRMatrix::laplacianMatrix(graph);
    Lamg<CSRMatrix> lamg;
    lamg.setupConnected(L);

    Vector rhs(201, 0.0);
    Vector result(201, 0.0);
    rhs[0] = -1.0;
    rhs[1] = 1.0;

    const auto status = lamg.solve(rhs, result);
    EXPECT_TRUE(vector_almost_equal(L * result, rhs))
        << "Lamg result: " << result << "\nL * result: " << L * result << "\nrhs: " << rhs;
    EXPECT_GT(result.length(), 0.1);
    EXPECT_LT(status.residual, 1e-6);
    EXPECT_TRUE(status.converged);
}

TEST(LAMGGTest, starEliminatedToOneNodeNonzeroGuess) {
    omp_set_num_threads(1);

    Graph graph(201);
    for (node v = 1; v < 201; ++v)
        graph.addEdge(0, v);

    const auto L = CSRMatrix::laplacianMatrix(graph);
    Lamg<CSRMatrix> lamg;
    lamg.setupConnected(L);

    Vector rhs(201, 0.0);
    rhs[0] = -1.0;
    rhs[1] = 1.0;

    // Nonzero guess: solveThread sets the reduction to tol*||b||/||L*x0-b||.
    // The one-node result is still the gauged solution, independent of x0.
    Vector result(201, 0.0);
    result[0] = 1e8;
    const double initialResidual = (rhs - L * result).length();
    EXPECT_GT(initialResidual, rhs.length() * 1e6);

    const auto status = lamg.solve(rhs, result);
    EXPECT_TRUE(vector_almost_equal(L * result, rhs))
        << "Lamg result: " << result << "\nL * result: " << L * result << "\nrhs: " << rhs;
    EXPECT_LT(status.residual, 1e-6 * rhs.length());
    EXPECT_TRUE(status.converged);
}

TEST(LAMGGTest, starEliminatedToOneNodeInexactRhsNonzeroGuess) {
    omp_set_num_threads(1);

    Graph graph(201);
    for (node v = 1; v < 201; ++v)
        graph.addEdge(0, v);

    const auto L = CSRMatrix::laplacianMatrix(graph);
    Lamg<CSRMatrix> lamg;
    lamg.setupConnected(L);

    // 0.1, 0.2 and 0.3 are not exact in binary, so the residual stays near 6e-17
    // instead of 0. With a zero residual, any threshold would report converged.
    Vector rhs(201, 0.0);
    rhs[0] = -0.3;
    rhs[5] = 0.1;
    rhs[7] = 0.2;

    // Far from the solution, a threshold scaled by ||b|| instead of ||L*x0-b||
    // drops below that residual.
    Vector result(201, 0.0);
    result[0] = 1e8;

    const auto status = lamg.solve(rhs, result);
    EXPECT_TRUE(vector_almost_equal(L * result, rhs))
        << "Lamg result: " << result << "\nL * result: " << L * result << "\nrhs: " << rhs;
    EXPECT_LT(status.residual, 1e-6 * rhs.length());
    EXPECT_TRUE(status.converged);
}

} /* namespace NetworKit */
