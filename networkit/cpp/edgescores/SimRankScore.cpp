/*  SimRankScore.cpp
 *
 *  Created on: 01.05.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/edgescores/SimRankScore.hpp>

namespace NetworKit {

SimRankScore::SimRankScore(const Graph &G, double similarityPropagationFactor, count maxIterations,
                           double tolerance)
    : EdgeScore<double>(G), similarityPropagationFactor(similarityPropagationFactor),
      maxIterations(maxIterations), tolerance(tolerance) {
    if (similarityPropagationFactor < 0.0 || similarityPropagationFactor > 1.0) {
        throw std::invalid_argument("similarityPropagationFactor must be in the range [0,1]");
    }

    if (maxIterations == 0) {
        throw std::invalid_argument("maxIterations must be greater than 0");
    }

    if (tolerance < 0.0) {
        throw std::invalid_argument("tolerance must be greater than or equal to 0");
    }
}

void SimRankScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("Graph must have edge IDs for SimRank calculation");
    }

    const index n = G->upperNodeIdBound();

    // This is practically unreachable in ordinary tests because constructing a graph
    // with n >= 2^32 nodes usually fails first, but the check prevents n * n from
    // wrapping before allocating the dense SimRank matrix.
    if (n != 0 && n > std::numeric_limits<index>::max() / n) {
        throw std::overflow_error("SimRankScore matrix size overflows");
    }

    const index matrixSize = n * n;

    auto matrixIndex = [n](node u, node v) -> index { return u * n + v; };

    std::vector<double> oldScore(matrixSize, 0.0);
    std::vector<double> newScore(matrixSize, 0.0);

    auto initDiagonal = [&](std::vector<double> &score) {
        G->parallelForNodes([&](node u) { score[matrixIndex(u, u)] = 1.0; });
    };

    initDiagonal(oldScore);

    double maxDifference = std::numeric_limits<double>::max();
    std::vector<double> threadMaxDifferences(omp_get_max_threads());

    for (count iterations = 0; iterations < maxIterations; ++iterations) {
        std::ranges::fill(newScore.begin(), newScore.end(), 0.0);
        std::ranges::fill(threadMaxDifferences, 0.0);

        initDiagonal(newScore);

        G->parallelForNodePairs([&](node u, node v) {
            const count degreeU = G->isDirected() ? G->degreeIn(u) : G->degree(u);
            const count degreeV = G->isDirected() ? G->degreeIn(v) : G->degree(v);

            double value = 0.0;

            if (degreeU > 0 && degreeV > 0) {
                double sum = 0.0;

                if (G->isDirected()) {
                    for (node a : G->inNeighborRange(u)) {
                        for (node b : G->inNeighborRange(v)) {
                            sum += oldScore[matrixIndex(a, b)];
                        }
                    }
                } else {
                    for (node a : G->neighborRange(u)) {
                        for (node b : G->neighborRange(v)) {
                            sum += oldScore[matrixIndex(a, b)];
                        }
                    }
                }

                value = similarityPropagationFactor * sum
                        / (static_cast<double>(degreeU) * static_cast<double>(degreeV));
            }

            const auto uv = matrixIndex(u, v);
            const auto vu = matrixIndex(v, u);

            newScore[uv] = value;
            newScore[vu] = value;

            const double difference =
                std::max(std::abs(value - oldScore[uv]), std::abs(value - oldScore[vu]));

            const int tid = omp_get_thread_num();
            threadMaxDifferences[tid] = std::max(threadMaxDifferences[tid], difference);
        });

        const double iterationMaxDifference = std::ranges::max(threadMaxDifferences);

        oldScore.swap(newScore);
        maxDifference = iterationMaxDifference;
        ++iterations;

        if (maxDifference < tolerance) {
            break;
        }
    }

    scoreData.assign(G->upperEdgeIdBound(), 0.0);

    G->parallelForEdges(
        [&](node u, node v, edgeid eid) { scoreData[eid] = oldScore[matrixIndex(u, v)]; });

    hasRun = true;
}

} // namespace NetworKit
