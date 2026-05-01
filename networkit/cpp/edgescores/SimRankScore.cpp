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
    double maxDifference = std::numeric_limits<double>::max();

    auto matrixIndex = [n](node u, node v) -> index { return u * n + v; };

    std::vector<double> oldScore(matrixSize, 0.0);
    std::vector<double> newScore(matrixSize, 0.0);

    G->forNodes([&](node u) { oldScore[matrixIndex(u, u)] = 1.0; });

    iterations = 0;
    for (count iter = 0; iter < maxIterations; ++iter) {
        std::fill(newScore.begin(), newScore.end(), 0.0);

        G->forNodes([&](node u) { newScore[matrixIndex(u, u)] = 1.0; });
        double iterationMaxDifference = 0.0;
        G->forNodes([&](node u) {
            const count degreeU = G->degreeIn(u);
            G->forNodes([&](node v) {
                if (u == v) {
                    return;
                }
                const count degreeV = G->degreeIn(v);
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
                    value =
                        similarityPropagationFactor * sum / static_cast<double>(degreeU * degreeV);
                }
                const auto index = matrixIndex(u, v);
                newScore[index] = value;
                iterationMaxDifference =
                    std::max(iterationMaxDifference, std::abs(value - oldScore[index]));
            });
        });
        oldScore.swap(newScore);
        maxDifference = iterationMaxDifference;
        iterations++;

        if (maxDifference < tolerance) {
            break;
        }
    }
    scoreData.assign(G->upperEdgeIdBound(), 0.0);
    G->forEdges([&](node u, node v, edgeid eid) { scoreData[eid] = oldScore[matrixIndex(u, v)]; });
    hasRun = true;
}

} // namespace NetworKit
