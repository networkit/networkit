/*
 * PageRank.cpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

// networkit-format

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/PageRank.hpp>

namespace NetworKit {

PageRank::PageRank(const Graph &G, double damp, double tol)
    : Centrality(G, true), damp(damp), tol(tol) {}

void PageRank::run() {
    Aux::SignalHandler handler;
    const auto n = G.numberOfNodes();
    const auto z = G.upperNodeIdBound();

    const auto teleportProb = (1.0 - damp) / static_cast<double>(n);
    scoreData.resize(z, 1.0 / static_cast<double>(n));
    std::vector<double> pr = scoreData;

    std::vector<double> deg(z, 0.0);
    G.parallelForNodes([&](const node u) { deg[u] = static_cast<double>(G.weightedDegree(u)); });

    auto converged([&]() {
        double diff = G.parallelSumForNodes([&](node u) {
            double d = scoreData[u] - pr[u];
            return d * d;
        });
        return (sqrt(diff) <= tol);
    });

    bool isConverged = false;
    while (!isConverged) {
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

        isConverged = converged();
        scoreData = pr;
    }

    handler.assureRunning();

    const auto sum = G.parallelSumForNodes([&](const node u) { return scoreData[u]; });

    // make sure scoreData sums up to 1
    assert(!Aux::NumericTools::equal(sum, 0.0, 1e-15));
    G.parallelForNodes([&](const node u) { scoreData[u] /= sum; });

    hasRun = true;
}

double PageRank::maximum() {
    return 1.0; // upper bound, could be tighter by assuming e.g. a star graph with n nodes
}

} /* namespace NetworKit */
