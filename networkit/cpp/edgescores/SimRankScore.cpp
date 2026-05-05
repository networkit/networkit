/*  SimRankScore.cpp
 *
 *  Created on: 01.05.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/auxiliary/NumericTools.hpp>
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

    const auto matrixIndex = [n](node u, node v) -> index { return u * n + v; };

    std::vector<double> oldScore(matrixSize, 0.0);
    std::vector<double> newScore(matrixSize, 0.0);

    const auto initDiagonal = [&](std::vector<double> &score) {
        G->parallelForNodes([&](node u) { score[matrixIndex(u, u)] = 1.0; });
    };

    initDiagonal(oldScore);

    std::vector<double> threadMaxDifferences(omp_get_max_threads());

    const auto computeSum = [&](node u, node v, auto &&getRange) -> double {
        double sum = 0.0;
        for (node a : getRange(u)) {
            for (node b : getRange(v)) {
                sum += oldScore[matrixIndex(a, b)];
            }
        }
        return sum;
    };

    const auto updateScoreAndDifference = [&](node u, node v, double sum, count degreeU,
                                              count degreeV) {
        double value = 0.0;
        if (degreeU > 0 && degreeV > 0) {
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
    };

    for (count iterations = 0; iterations < maxIterations; ++iterations) {
        std::ranges::fill(newScore, 0.0);
        std::ranges::fill(threadMaxDifferences, 0.0);

        initDiagonal(newScore);
        if (G->isDirected()) {
            G->parallelForNodePairs([&](node u, node v) {
                const double sum = computeSum(u, v, [&](node w) { return G->inNeighborRange(w); });
                updateScoreAndDifference(u, v, sum, G->degreeIn(u), G->degreeIn(v));
            });
        } else {
            G->parallelForNodePairs([&](node u, node v) {
                const double sum = computeSum(u, v, [&](node w) { return G->neighborRange(w); });
                updateScoreAndDifference(u, v, sum, G->degree(u), G->degree(v));
            });
        }

        oldScore.swap(newScore);

        const double iterationMaxDifference = std::ranges::max(threadMaxDifferences);
        if (iterationMaxDifference < tolerance) {
            break;
        }
    }

    scoreData.assign(G->upperEdgeIdBound(), 0.0);

    G->parallelForEdges(
        [&](node u, node v, edgeid eid) { scoreData[eid] = oldScore[matrixIndex(u, v)]; });

    hasRun = true;
}

} // namespace NetworKit
