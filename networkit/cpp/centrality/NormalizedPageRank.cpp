/*
 * PageRank.cpp
 *
 *  Created on: 03.03.2021
 *      Author: Petersen
 */

// networkit-format

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/NormalizedPageRank.hpp>

namespace NetworKit {

NormalizedPageRank::NormalizedPageRank(const Graph &G, double damp, double tol)
    : Centrality(G, true), damp(damp), tol(tol) {}

void NormalizedPageRank::run() {
    Aux::SignalHandler handler;
    const auto n = G.numberOfNodes();
    const auto z = G.upperNodeIdBound();

    const auto teleportProb = (1.0 - damp) / static_cast<double>(n);
    scoreData.resize(z, 1.0 / static_cast<double>(n));
    std::vector<double> pr = scoreData;

    std::vector<double> deg(z, 0.0);
    G.parallelForNodes([&](const node u) { deg[u] = static_cast<double>(G.weightedDegree(u)); });

    iterations = 0;

    auto sumL1Norm = [&](const node u) { return std::abs(scoreData[u] - pr[u]); };

    auto sumL2Norm = [&](const node u) {
        const auto d = scoreData[u] - pr[u];
        return d * d;
    };

    auto converged([&]() {
        if (iterations >= maxIterations) {
            return true;
        }

        if (norm == Norm2::L2Norm) {
            return std::sqrt(G.parallelSumForNodes(sumL2Norm)) <= tol;
        }

        return G.parallelSumForNodes(sumL1Norm) <= tol;
    });

    bool isConverged = false;
    do {
        handler.assureRunning();
        G.balancedParallelForNodes([&](const node u) {
            pr[u] = 0.0;
            G.forInEdgesOf(u, [&](const node u, const node v, const edgeweight w) {
                // note: inconsistency in definition in Newman's book (Ch. 7) regarding directed
                // graphs we follow the verbal description, which requires to sum over the incoming
                // edges
                pr[u] += scoreData[v] * w / deg[v];
            });
            pr[u] *= damp;
            pr[u] += teleportProb;
        });

        ++iterations;
        isConverged = converged();
        std::swap(pr, scoreData);
    } while (!isConverged);

    handler.assureRunning();

    // calculate sum of the Pagerank of dangling Nodes for normalization
    const auto sum = G.parallelSumForNodes([&](const node u) {
        if (G.degree(u) == 0.0)
            return scoreData[u];
        return 0.0;
    });

    const auto normFactor = (1.0 / n) * (tol + ((1.0 - tol) * sum));

    G.parallelForNodes([&](const node u) { scoreData[u] /= normFactor; });

    hasRun = true;
}

// Unlike regular Page Rank there isn't a universal upper bound. So we need to iterate.
double NormalizedPageRank::maximum() {
    auto max = scoreData[0];
    G.balancedParallelForNodes([&](const node u) {
        if (max < scoreData[u])
            max = scoreData[u];
    });
    return max;
}

} /* namespace NetworKit */
