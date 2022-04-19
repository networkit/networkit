#include <omp.h>
#include <networkit/centrality/LocalSquareClusteringCoefficient.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

namespace NetworKit {

LocalSquareClusteringCoefficient::LocalSquareClusteringCoefficient(const Graph &G)
    : Centrality(G, false, false) {
    if (G.isDirected())
        throw std::runtime_error("Not implemented: Local square clustering coefficient is "
                                 "currently not implemented for directed graphs");
    if (G.numberOfSelfLoops())
        throw std::runtime_error(
            "Local square clustering coefficient implementation does not support graphs with "
            "self-loops. Call Graph.removeSelfLoops() first.");
}

void LocalSquareClusteringCoefficient::run() {
    count z = G.upperNodeIdBound();
    scoreData.clear();
    scoreData.resize(z);

    G.balancedParallelForNodes([&](const node u) {
        double squares = 0;           // Number of squares found.
        double potential_squares = 0; // Maximum number of possible squares.
        const auto neighborsU = G.neighborRange(u);

        // Iterate over all combinations of neighbors.
        for (auto iterV = neighborsU.begin(); iterV != neighborsU.end(); iterV++) {
            for (auto iterW = std::next(iterV); iterW != neighborsU.end(); ++iterW) {
                // Find the number of common neighbors (including the node `u`) and count squares.
                const auto numCommonNeighbors =
                    NeighborhoodUtility::getCommonNeighbors(G, *iterV, *iterW).size();
                squares += numCommonNeighbors - 1;
                // Update the number of squares that could be realized.
                potential_squares += G.degree(*iterV) + G.degree(*iterW) - numCommonNeighbors - 1;
                if (G.hasEdge(*iterV, *iterW)) {
                    potential_squares -= 2;
                }
            }
        }
        if (potential_squares > 0) {
            squares /= potential_squares;
        }
        scoreData[u] = squares;
    });

    hasRun = true;
}

double LocalSquareClusteringCoefficient::maximum() {
    return 1.0;
}

} // namespace NetworKit
