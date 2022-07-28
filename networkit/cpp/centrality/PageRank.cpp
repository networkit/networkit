/*
 * PageRank.cpp
 *
 *  Created on: 19.03.2014
 *      Authors: Henning Meyerhenke
 *               Fabian Brandt-Tumescheit <brandtfa@hu-berlin.de>
 */

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/PageRank.hpp>

namespace NetworKit {

PageRank::PageRank(const Graph &G, double damp, double tol, bool normalized,
                   SinkHandling distributeSinks)
    : Centrality(G, true), damp(damp), tol(tol), normalized(normalized),
      distributeSinks(distributeSinks) {}

void PageRank::run() {
    Aux::SignalHandler handler;
    const auto n = G.numberOfNodes();
    const auto z = G.upperNodeIdBound();

    const auto teleportProb = (1.0 - damp) / static_cast<double>(n);
    const double factor = damp / static_cast<double>(n);
    scoreData.resize(z, 1.0 / static_cast<double>(n));
    std::vector<double> pr = scoreData;

    std::vector<double> deg(z, 0.0);
    G.parallelForNodes([&](const node u) { deg[u] = static_cast<double>(G.weightedDegree(u)); });

    std::vector<node> sinks;
    if (G.isDirected() && ((distributeSinks == SinkHandling::DISTRIBUTE_SINKS) || normalized)) {
        G.forNodes([&](const node u) {
            if (G.degree(u) == 0) {
                sinks.push_back(u);
            }
        });
    }
    count nSinks = sinks.size();

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

        if (norm == Norm::L2_NORM) {
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

        // For directed graphs sink-handling is needed to fulfill |pr| == 1 in each step. Otherwise
        // probability mass would be leaked, creating wrong results. For this, we add edges from
        // sinks to all other nodes. This is described amongst others in "PageRank revisited."
        // by M. Brinkmeyer et al. (2005).
        if (G.isDirected() && ((distributeSinks == SinkHandling::DISTRIBUTE_SINKS) || normalized)) {
            double totalSinkContrib = 0.0;
#pragma omp parallel for reduction(+ : totalSinkContrib)
            for (omp_index i = 0; i < static_cast<omp_index>(nSinks); i++) {
                totalSinkContrib += factor * scoreData[sinks[i]];
            }
            G.balancedParallelForNodes([&](const node u) { pr[u] += totalSinkContrib; });
        }

        ++iterations;
        isConverged = converged();
        std::swap(pr, scoreData);
    } while (!isConverged);

    handler.assureRunning();

    // Post-processing for normalized PageRank
    if (normalized) {
        double normFactor;
        if (G.isDirected()) {
            // Calculate sum of dangling Nodes for normalization
            double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
            for (omp_index i = 0; i < static_cast<omp_index>(nSinks); i++) {
                sum += scoreData[sinks[i]];
            }
            normFactor = (1.0 / static_cast<double>(n)) * ((1.0 - damp) + (damp * sum));
        } else {
            normFactor = teleportProb;
        }
        G.parallelForNodes([&](const node u) { scoreData[u] /= normFactor; });

        // Post-processing for non-normalized PageRank
    } else {
        if (G.isDirected() && distributeSinks == SinkHandling::NO_SINK_HANDLING) {
            // In case no sink handling was done, make sure that |pr| == 1
            const auto sum = G.parallelSumForNodes([&](const node u) { return scoreData[u]; });
            G.parallelForNodes([&](const node u) { scoreData[u] /= sum; });
        }
    }
    // calculate the maxium
    max = scoreData[0];
    G.balancedParallelForNodes([&](const node u) { Aux::Parallel::atomic_max(max, scoreData[u]); });
    hasRun = true;
}

double PageRank::maximum() {
    return max.load(std::memory_order_relaxed);
}

} /* namespace NetworKit */
