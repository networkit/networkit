/*
 * Lamg.hpp
 *
 *  Created on: Oct 20, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_NUMERICS_LAMG_LAMG_HPP_
#define NETWORKIT_NUMERICS_LAMG_LAMG_HPP_

#include <omp.h>
#include <vector>

#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/components/ParallelConnectedComponents.hpp>
#include <networkit/numerics/GaussSeidelRelaxation.hpp>
#include <networkit/numerics/LAMG/LevelHierarchy.hpp>
#include <networkit/numerics/LAMG/MultiLevelSetup.hpp>
#include <networkit/numerics/LAMG/SolverLamg.hpp>
#include <networkit/numerics/LinearSolver.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 * Represents the interface to the Lean Algebraic Multigrid (LAMG) graph Laplacian linear solver
 * by Oren E. Livne and Achi Brandt.
 * @see Livne, Oren E., and Achi Brandt. "Lean algebraic multigrid (LAMG): Fast graph Laplacian
 * linear solver." SIAM Journal on Scientific Computing 34.4 (2012): B499-B522.
 */
template <class Matrix>
class Lamg : public LinearSolver<Matrix> {
private:
    // for a graph/matrix with K components and T threads:

    bool validSetup;
    const GaussSeidelRelaxation<Matrix> smoother;
    const MultiLevelSetup<Matrix> lamgSetup;
    Matrix laplacianMatrix;
    count numComponents;

    // maps a node id from global index (in L) to its node id in the component (block of L) such
    // that in each component, nodes are indexed from 0 to size(component)-1. Size n
    // not used in case there is only a single component
    std::vector<index> graph2Components;

    // for each component, a vector of ids that are in this component. Size K
    std::vector<std::vector<index>> components;
    // one LevelHierarchy for each component. Size K
    std::vector<LevelHierarchy<Matrix>> compHierarchies;

    // these are internal data structures used in the solve function
    // one object for each component and each thread - access: [thread][component]
    mutable std::vector<std::vector<SolverLamg<Matrix>>> compSolvers;
    // one object for each component and each thread - access: [thread][component]
    mutable std::vector<std::vector<LAMGSolverStatus>> compStati;

    // one object for each component and each thread - access: [thread][component]
    mutable std::vector<std::vector<Vector>> initialVectors;
    // one object for each component and each thread - access: [thread][component]
    mutable std::vector<std::vector<Vector>> rhsVectors;

    // initializing step 1
    // initializes numComponents, graph2components, components, compHierarchies
    // precondition: laplacianMatrix is set correctly
    void initializeOneComponent();
    void initializeMultipleComponents(const Graph &G, const ComponentDecomposition &decomp);

    // initializing step 2
    // initializes compSolvers, compStati, initialVectors, rhsVectors
    // Precondition: the following members are set correctly:
    // numComponents, components, compHierarchies
    void initializeInternalDatastructures() const;

    // solves on the given thread
    SolverStatus solveThread(const Vector &rhs, Vector &result, count maxConvergenceTime,
                             count maxIterations, index threadId) const;

public:
    /**
     * Construct a solver with the given @a tolerance. The relative residual ||Ax-b||/||b|| will be
     * less than or equal to
     * @a tolerance after the solver finished.
     * @param tolerance
     */
    Lamg(const double tolerance = 1e-6)
        : LinearSolver<Matrix>(tolerance), validSetup(false), lamgSetup(smoother),
          numComponents(0) {}
    /** Default destructor */
    ~Lamg() override = default;

    /**
     * Compute the multigrid hierarchy for the given Laplacian matrix @a laplacianMatrix.
     * @param laplacianMatrix
     * @note This method also works for disconnected graphs. If you know that the graph is
     * connected, it is faster to use @ref setupConnected instead.
     */
    void setup(const Matrix &laplacianMatrix) override;

    /**
     * Compute the multigrid hierarchy for the Laplacian matrix that
     * corresponds to the graph @a G.
     * @param laplacianMatrix
     * @param G Graph that corresponds to the @a laplacianMatrix
     */
    void setup(const Graph &graph) override;

    /**
     * Compute the multigrid hierarchy for the given Laplacian matrix @a laplacianMatrix that
     * corresponds to the graph @a G.
     * @param laplacianMatrix
     * @param G Graph that corresponds to the @a laplacianMatrix
     * @note This setup method allows you to skip some computation required for setting up LAMG. You
     * can provide your matrix and graph, which are required for setting up LAMG, via this method to
     * prevent duplicate computation. Note that the output is undefined if the two parameters do not
     * correspond the the same graph.
     */
    void setup(const Matrix &laplacianMatrix, const Graph &G);

    /**
     * Compute the multigrid hierarchy for the Laplacian matrix that
     * corresponds to the graph @a G using the existing component decomposition @a decomp.
     * @param G Graph that corresponds to the @a laplacianMatrix
     * @param decomp ComponentDecomposition for the graph @a G
     * @note This setup method allows you to skip some computation required for setting up LAMG. You
     * can provide your graph and decomposition, which are required for setting up LAMG, via this
     * method to prevent duplicate computation. Note that the output is undefined if the two
     * parameters do not correspond the the same graph.
     */
    void setup(const Graph &G, const ComponentDecomposition &decomp);

    /**
     * Compute the multigrid hierarchy for the given Laplacian matrix @a laplacianMatrix that
     * corresponds to the graph @a G using the existing component decomposition @a decomp.
     * @param laplacianMatrix
     * @param G Graph that corresponds to the @a laplacianMatrix
     * @param decomp ComponentDecomposition for the graph @a G
     * @note This setup method allows you to skip some computation required for setting up LAMG. You
     * can provide your graph and decomposition, which are required for setting up LAMG, via
     * this method to prevent duplicate computation. Note that the output is undefined if the three
     * parameters do not correspond the the same graph.
     */
    void setup(const Matrix &laplacianMatrix, const Graph &G, const ComponentDecomposition &decomp);

    /**
     * Compute the multigrid hierarchy for the given Laplacian matrix @a laplacianMatrix.
     * @param laplacianMatrix
     * @note The graph has to be connected for this method to work. Otherwise the output is
     * undefined.
     */
    void setupConnected(const Matrix &laplacianMatrix) override;

    /**
     * Computes the @a result for the matrix currently setup and the right-hand side @a rhs.
     * The maximum spent time can be specified by @a maxConvergenceTime and the maximum number of
     * iterations can be set by @a maxIterations.
     * @param rhs
     * @param result
     * @param maxConvergenceTime
     * @param maxIterations
     * @return A @ref SolverStatus object which provides some statistics like the final absolute
     * residual.
     */
    SolverStatus solve(const Vector &rhs, Vector &result, count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max()) const override;

    /**
     * Compute the @a results for the matrix currently setup and the right-hand sides @a rhs.
     * The maximum spent time for each system can be specified by @a maxConvergenceTime and the
     * maximum number of iterations can be set by @a maxIterations.
     * @param rhs
     * @param results
     * @param maxConvergenceTime
     * @param maxIterations
     */
    std::vector<SolverStatus>
    parallelSolve(const std::vector<Vector> &rhs, std::vector<Vector> &results,
                  count maxConvergenceTime = 5 * 60 * 1000,
                  count maxIterations = std::numeric_limits<count>::max()) const override;

    /**
     * Abstract parallel solve function that computes and processes results using @a resultProcessor
     * for the matrix currently setup and the right-hand sides (size of @a rhsSize) provided by @a
     * rhsLoader. The maximum spent time for each system can be specified by @a maxConvergenceTime
     * and the maximum number of iterations can be set by @a maxIterations.
     * @param rhsLoader
     * @param resultProcessor
     * @param rhsSize
     * @param maxConvergenceTime
     * @param maxIterations
     * @note If the solver does not support parallelism during solves, this function falls back to
     * solving the systems sequentially.
     * @note The @a rhsLoader and @a resultProcessor are called from multiple threads in parallel.
     */
    template <typename RHSLoader, typename ResultProcessor>
    void parallelSolve(const RHSLoader &rhsLoader, const ResultProcessor &resultProcessor,
                       std::pair<count, count> rhsSize, count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max()) const {

        const count n = rhsSize.first;
        const count m = rhsSize.second;
        const index numThreads = omp_get_max_threads();

        std::vector<Vector> results(numThreads, Vector(m));
        std::vector<Vector> RHSs(numThreads, Vector(m));

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
            const index threadId = omp_get_thread_num();

            const Vector &rhs = rhsLoader(i, RHSs[threadId]);
            Vector &result = results[threadId];

            solveThread(rhs, result, maxConvergenceTime, maxIterations, threadId);
            resultProcessor(i, result);
        }
    }
};

template <class Matrix>
void Lamg<Matrix>::initializeInternalDatastructures() const {
    // preconditions:
    // - numComponents is correct
    // - components vector is correct
    // - compHierarchies is setup

    // first, resize outer vectors
    // then iterate elements and initialize them correctly

    initialVectors.resize(omp_get_max_threads());
    rhsVectors.resize(omp_get_max_threads());
    compSolvers.resize(omp_get_max_threads());
    compStati.resize(omp_get_max_threads());

#pragma omp parallel for
    for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
        // resize inner vectors
        initialVectors[thread].resize(numComponents);
        rhsVectors[thread].resize(numComponents);
        compStati[thread].resize(numComponents);

        // there is no default constructor for SolverLamg - reserve memory and construct later
        compSolvers[thread].clear();
        compSolvers[thread].reserve(numComponents);

        // iterate components
        for (index compIdx = 0; compIdx < numComponents; ++compIdx) {
            auto &component = components[compIdx];

            // Vector objects do not have a resize method - create new objects instead
            initialVectors[thread][compIdx] = Vector(component.size());
            rhsVectors[thread][compIdx] = Vector(component.size());

            compSolvers[thread].emplace_back(compHierarchies[compIdx], smoother);

            LAMGSolverStatus status;
            status.desiredResidualReduction =
                this->tolerance * component.size() / laplacianMatrix.numberOfColumns();
            compStati[thread][compIdx] = status;
        }
    }
}

template <class Matrix>
void Lamg<Matrix>::initializeOneComponent() {
    numComponents = 1;
    components.resize(1);
    components[0].resize(laplacianMatrix.numberOfColumns());
    std::iota(components[0].begin(), components[0].end(), 0);
    graph2Components.clear();
    compHierarchies = std::vector<LevelHierarchy<Matrix>>(1);
    lamgSetup.setup(laplacianMatrix, compHierarchies[0]);
    initializeInternalDatastructures();
}

template <class Matrix>
void Lamg<Matrix>::initializeMultipleComponents(const Graph &G,
                                                const ComponentDecomposition &decomp) {

    numComponents = decomp.numberOfComponents();

    components.resize(numComponents);
    graph2Components.resize(G.numberOfNodes());

    compHierarchies.resize(numComponents);

    index compIdx = 0;
    for (const auto &partitionComponent : decomp.getPartition().getSubsets()) {
        components[compIdx].resize(partitionComponent.size());
        std::copy(partitionComponent.begin(), partitionComponent.end(),
                  components[compIdx].begin());

        auto &component = components[compIdx];

        index idx = 0;
        for (node u : component) {
            graph2Components[u] = idx;
            idx++;
        }

        // it is not clear to me if LevelHierarchy can be reused for different matrices. Hence,
        // construct new objects
        compHierarchies[compIdx] = LevelHierarchy<Matrix>();

        // construct block matrix via neighborsOf in G
        std::vector<Triplet> triplets;
        for (node u : component) {
            G.forNeighborsOf(u, [&](node v, edgeweight w) {
                triplets.push_back({graph2Components[u], graph2Components[v], w});
            });
        }
        Matrix compMatrix(component.size(), component.size(), triplets);
        lamgSetup.setup(compMatrix, compHierarchies[compIdx]);

        ++compIdx;
    }
    initializeInternalDatastructures();
}

template <class Matrix>
void Lamg<Matrix>::setupConnected(const Matrix &laplacianMatrix) {
    assert([&]() {
        auto G = MatrixTools::matrixToGraph(laplacianMatrix);
        ParallelConnectedComponents cc(G);
        cc.run();
        return cc.numberOfComponents() == 1;
    }());
    this->laplacianMatrix = laplacianMatrix;
    initializeOneComponent();
    validSetup = true;
}

template <class Matrix>
void Lamg<Matrix>::setup(const Matrix &laplacianMatrix, const Graph &G,
                         const ComponentDecomposition &decomp) {
    this->laplacianMatrix = laplacianMatrix;
    numComponents = decomp.numberOfComponents();
    if (numComponents == 1) {
        initializeOneComponent();
    } else {
        initializeMultipleComponents(G, decomp);
    }
    validSetup = true;
}

template <class Matrix>
void Lamg<Matrix>::setup(const Graph &G, const ComponentDecomposition &decomp) {
    setup(Matrix::laplacianMatrix(G), G, decomp);
}

template <class Matrix>
void Lamg<Matrix>::setup(const Graph &G) {
    setup(Matrix::laplacianMatrix(G), G);
}

template <class Matrix>
void Lamg<Matrix>::setup(const Matrix &laplacianMatrix, const Graph &G) {
    ParallelConnectedComponents con(G, false);
    con.run();
    setup(laplacianMatrix, G, con);
}

template <class Matrix>
void Lamg<Matrix>::setup(const Matrix &laplacianMatrix) {
    Graph G = MatrixTools::matrixToGraph(laplacianMatrix);
    setup(laplacianMatrix, G);
}

template <class Matrix>
SolverStatus Lamg<Matrix>::solveThread(const Vector &rhs, Vector &result, count maxConvergenceTime,
                                       count maxIterations, const index threadId) const {
    if (!validSetup)
        throw std::runtime_error("LAMG is not properly setup!");
    if (result.getDimension() != laplacianMatrix.numberOfColumns()
        || rhs.getDimension() != laplacianMatrix.numberOfRows()) {
        throw std::runtime_error("Wrong matrix dimensions for given vectors.");
    }
    SolverStatus status;

    if (numComponents == 1) {
        LAMGSolverStatus &stat = compStati[threadId][0];
        stat.desiredResidualReduction =
            this->tolerance * rhs.length() / (laplacianMatrix * result - rhs).length();
        stat.maxIters = maxIterations;
        stat.maxConvergenceTime = maxConvergenceTime;
        compSolvers[threadId][0].solve(result, rhs, stat);

        status.residual = stat.residual;
        status.numIters = stat.numIters;
        status.converged = stat.converged;
    } else {
        // solve on every component
        count maxIters = 0;
        for (index componentId = 0; componentId < components.size(); ++componentId) {
            LAMGSolverStatus &stat = compStati[threadId][componentId];
            Vector &componentRhs = rhsVectors[threadId][componentId];
            Vector &componentResult = initialVectors[threadId][componentId];

            // setup component vectors
            for (auto element : components[componentId]) {
                componentResult[graph2Components[element]] = result[element];
                componentRhs[graph2Components[element]] = rhs[element];
            }

            // setup solver status
            double resReduction =
                this->tolerance * componentRhs.length()
                / (compHierarchies[componentId].at(0).getLaplacian() * componentResult
                   - componentRhs)
                      .length();
            stat.desiredResidualReduction =
                resReduction * components[componentId].size() / laplacianMatrix.numberOfRows();
            stat.maxIters = maxIterations;
            stat.maxConvergenceTime = maxConvergenceTime;

            compSolvers[threadId][componentId].solve(componentResult, componentRhs, stat);
            // write solution back to result
            for (auto element : components[componentId]) {
                result[element] = componentResult[graph2Components[element]];
            }

            maxIters = std::max(maxIters, stat.numIters);
        }

        status.residual = (rhs - laplacianMatrix * result).length();
        status.converged = status.residual <= this->tolerance;
        status.numIters = maxIters;
    }

    return status;
}

template <class Matrix>
SolverStatus Lamg<Matrix>::solve(const Vector &rhs, Vector &result, count maxConvergenceTime,
                                 count maxIterations) const {
    return solveThread(rhs, result, maxConvergenceTime, maxIterations, 0);
}

template <class Matrix>
std::vector<SolverStatus>
Lamg<Matrix>::parallelSolve(const std::vector<Vector> &rhs, std::vector<Vector> &results,
                            count maxConvergenceTime, count maxIterations) const {
    std::vector<SolverStatus> stati(rhs.size());

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(rhs.size()); ++i) {
        const index threadId = omp_get_thread_num();
        stati[i] = solveThread(rhs[i], results[i], maxConvergenceTime, maxIterations, threadId);
    }

    return stati;
}

extern template class Lamg<CSRMatrix>;
extern template class Lamg<DenseMatrix>;
extern template class Lamg<DynamicMatrix>;

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_LAMG_HPP_
