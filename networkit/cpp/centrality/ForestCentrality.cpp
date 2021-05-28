/*
 * ForestCentrality.cpp
 *
 *  Created on: 17.10.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/ForestCentrality.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/Preconditioner/DiagonalPreconditioner.hpp>

namespace NetworKit {

ForestCentrality::ForestCentrality(const Graph &G, node root, double inputEpsilon, double kappa)
    : Centrality(G), root(root), epsilon(inputEpsilon), kappa(kappa), volG(GraphTools::volume(G)),
      numberOfUSTs(computeNumberOfUSTs()), decDegree(G.nodeRange().begin(), G.nodeRange().end()) {

    if (G.isWeighted())
        WARN("ForestCentrality ignores edge weights.");

    if (G.isDirected())
        throw std::runtime_error("The input graph must be undirected.");

    if (G.numberOfNodes() != G.upperNodeIdBound())
        throw std::runtime_error("The graph must be compact; use GraphTools::getCompactedGraph");

    if (G.degree(root) != G.numberOfNodes() - 1)
        throw std::runtime_error("The input graph is not an augmented graph. Create an augmented "
                                 "graph with GraphTools::createAugmentedGraph.");

    const auto n = G.upperNodeIdBound();
    parentsGlobal.resize(omp_get_max_threads(), std::vector<node>(n));
    statusGlobal.resize(omp_get_max_threads(), std::vector<uint8_t>(n));
    std::sort(decDegree.begin(), decDegree.end(),
              [&G](node n1, node n2) { return G.degree(n1) < G.degree(n2); });

    diagonal.resize(n);

    approxEffResistanceGlobal.resize(omp_get_max_threads(), std::vector<count>(n));
    G.forNodes([&](node u) {
        const auto degU = G.degree(u);
        if (degU == 0)
            uniformDistr.emplace_back();
        else
            uniformDistr.emplace_back(0, degU - 1);
    });
}

void ForestCentrality::run() {
    sampleUSTs();
    computeDiagonal();
    computeScores();
    hasRun = true;
}

void ForestCentrality::sampleUSTs() {
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0)
            solveLinearSystem(); // Thread 0 solves the linear system

#pragma omp for schedule(dynamic)
        for (omp_index i = 0; i < static_cast<omp_index>(numberOfUSTs); ++i) {
            // Getting thread-local vectors
            const omp_index thread = omp_get_thread_num();
            auto &parent = parentsGlobal[thread];
            auto &r = approxEffResistanceGlobal[thread];
            auto &status = statusGlobal[thread];
            auto &generator = Aux::Random::getURNG();
            uint8_t rootStatus = ++status[root];

            count nodesInUST = 1;

            for (const node walkStart : decDegree) {
                if (status[walkStart] == rootStatus)
                    continue;

                node currentNode = walkStart;
                do {
                    auto idx = uniformDistr[currentNode](generator);
                    assert(idx < G.degree(currentNode));
                    const node randomNeighbor = G.getIthNeighbor(currentNode, idx);
                    parent[currentNode] = randomNeighbor;
                    currentNode = randomNeighbor;
                } while (status[currentNode] != rootStatus);

                const node walkEnd = currentNode;
                node next = parent[walkStart];
                currentNode = walkStart;
                do {
                    status[currentNode] = rootStatus;
                    r[currentNode] += static_cast<count>(next == root);
                    currentNode = next;
                    next = parent[next];
                    ++nodesInUST;
                } while (currentNode != walkEnd);

                if (nodesInUST == G.numberOfNodes())
                    break;
            }
        }
    }
}

void ForestCentrality::solveLinearSystem() {
    const auto n = static_cast<double>(G.numberOfNodes());
    const auto L = CSRMatrix::laplacianMatrix(G);
    const double c = n / volG;
    const double eta = kappa * epsilon / (6. * (c + 2.) * volG);
    ConjugateGradient<CSRMatrix, DiagonalPreconditioner> cg(eta);
    cg.setupConnected(L);

    Vector rhs(G.numberOfNodes());
    linearSysSol = Vector(G.numberOfNodes());
    G.parallelForNodes([&](node u) { rhs[u] = -1. / n; });
    rhs[root] += 1.;
    cg.solve(rhs, linearSysSol);

    // Ensure that column sum is 0
    const double columnSum = G.parallelSumForNodes([&](node u) { return linearSysSol[u]; });
    linearSysSol -= columnSum / n;
}

void ForestCentrality::computeDiagonal() {
    diagonal.resize(G.upperNodeIdBound());
    G.parallelForNodes([&](node u) {
        double aggr = 0;
        for (const auto &approxEffResLocal : approxEffResistanceGlobal)
            aggr += static_cast<double>(approxEffResLocal[u]);
        diagonal[u] = std::max(0., aggr / static_cast<double>(numberOfUSTs) - linearSysSol[root]
                                       + 2. * linearSysSol[u]);
    });
}

void ForestCentrality::computeScores() {
    const auto n = static_cast<double>(G.numberOfNodes());
    scoreData.resize(G.numberOfNodes());
    const double trace = std::accumulate(diagonal.begin(), diagonal.end(), 0.0);
    G.parallelForNodes([&](node u) {
        const double farness = n * diagonal[u] + trace - 2.0;
        scoreData[u] = n / farness;
    });
}

} // namespace NetworKit
