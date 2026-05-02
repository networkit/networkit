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
      maxIterations(maxIterations), tolerance(tolerance), iterations(0) {
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

    if (n != 0 && n > std::numeric_limits<index>::max() / n) {
        throw std::overflow_error("SimRankScore matrix size overflows");
    }

    const index matrixSize = n * n;

    auto matrixIndex = [n](node u, node v) -> index {
        return u * n + v;
    };

    std::vector<double> oldScore(matrixSize, 0.0);
    std::vector<double> newScore(matrixSize, 0.0);

    G->parallelForNodes([&](node u) {
        oldScore[matrixIndex(u, u)] = 1.0;
    });

    iterations = 0;
    double maxDifference = std::numeric_limits<double>::max();

    for (count iter = 0; iter < maxIterations; ++iter) {
        std::fill(newScore.begin(), newScore.end(), 0.0);

        G->parallelForNodes([&](node u) {
            newScore[matrixIndex(u, u)] = 1.0;
        });

        std::vector<double> threadMaxDifferences(omp_get_max_threads(), 0.0);

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

            const double difference = std::max(
                std::abs(value - oldScore[uv]),
                std::abs(value - oldScore[vu])
            );

            const int tid = omp_get_thread_num();
            threadMaxDifferences[tid] =
                std::max(threadMaxDifferences[tid], difference);
        });

        double iterationMaxDifference = 0.0;
        for (const double localDifference : threadMaxDifferences) {
            iterationMaxDifference =
                std::max(iterationMaxDifference, localDifference);
        }

        oldScore.swap(newScore);
        maxDifference = iterationMaxDifference;
        ++iterations;

        if (maxDifference < tolerance) {
            break;
        }
    }

    scoreData.assign(G->upperEdgeIdBound(), 0.0);

    G->parallelForEdges([&](node u, node v, edgeid eid) {
        scoreData[eid] = oldScore[matrixIndex(u, v)];
    });

    hasRun = true;
}

} // namespace NetworKit
