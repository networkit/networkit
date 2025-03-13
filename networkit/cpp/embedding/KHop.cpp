/*
 *  KHop.cpp
 *
 *
 *  Created on: 23.06.2024
 *      Author: Alexander K. Ziebs
 *
 */

#include <tuple>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/linkprediction/OverlapCoefficient.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"
#include "LearnEmbedding.hpp"
#include <networkit/embedding/KHop.hpp>

namespace NetworKit {

KHop::KHop(const Graph &G, size_t K, count L, count N, count D, KHopMode M, count winSize,
           count iterations)
    : G(G), K(K), L(L), N(N), D(D), M(M), winSize(winSize), iterations(iterations) {

    if (L < 1)
        throw std::runtime_error("L must be greater than 0.");
    if (N < 1)
        throw std::runtime_error("N must be greater than 0.");
    if (D < 1)
        throw std::runtime_error("D must be greater than 0.");
    if (K < 1)
        throw std::runtime_error("K must be greater than 0.");
    if (winSize < 1)
        throw std::runtime_error("winSize must be greater than 0.");
    if (iterations < 1)
        throw std::runtime_error("iterations must be greater than 0.");

    if (G.numberOfNodes() != G.upperNodeIdBound()) {
        throw std::runtime_error("The node ids of the graph must be continuous.");
    }

    bool hasIsolatedNodes = false;
    for (node u : G.nodeRange()) {
        if (G.isIsolated(u)) {
            hasIsolatedNodes = true;
            break;
        }
    }
    if (hasIsolatedNodes)
        throw std::runtime_error("Isolated nodes are not allowed.");

    // Initialize k-hop graph
    computeKHopGraph();
    generateEdgeWeights();
}

void KHop::run() {
    using namespace Embedding;

    Aux::SignalHandler handler;

    TRACE("preprocess transition probabilities ...");

    BiasedRandomWalk brw(&kHopGraph);
    brw.preprocessTransitionProbs(1, 1);
    handler.assureRunning();

    TRACE("do biased walks ...");
    auto walks = brw.doWalks(L, N);
    handler.assureRunning();

    TRACE("learn embeddings ...");
    features = learnEmbeddings(walks, kHopGraph.numberOfNodes(), D, winSize, iterations);
    handler.assureRunning();

    hasRun = true;
}

const std::vector<std::vector<float>> &KHop::getFeatures() const {
    assureFinished();
    return features;
}

void KHop::generateEdgeWeights() {
    OverlapCoefficient oCoeff;
    oCoeff.setGraph(kHopGraph);

    kHopGraph.parallelForNodes([&](node u) {
        count neighborIndex = 0;
        kHopGraph.forNeighborsOf(u, [&](node v) {
            neighborIndex++;
            if (v > u) {
                kHopGraph.setWeightAtIthNeighbor(unsafe, u, neighborIndex,
                                                 std::exp(oCoeff.run(u, v)));
            }
        });
    });
}

void KHop::computeKHopGraph() {
    kHopGraph = NetworKit::Graph(G.numberOfNodes(), true, G.isDirected());

    auto prunedBFS = [&](node startNode, count bfsLevel, bool reverse) {
        std::vector<bool> visited(G.numberOfNodes(), 0);
        std::fill(visited.begin(), visited.end(), false);
        visited[startNode] = true;
        count currentLevel = 0;

        std::queue<node> q0, q1;
        q0.push(startNode);

        auto processNeighbor = [&](node v) -> void {
            if (visited[v])
                return;
            visited[v] = true;
            q1.push(v);
            if (M == KHopMode::STRICT && currentLevel < bfsLevel)
                return;
            if (reverse)
                kHopGraph.addPartialInEdge(unsafe, startNode, v);
            else
                kHopGraph.addPartialOutEdge(unsafe, startNode, v);
        };

        do {
            currentLevel++;
            do {
                const node u = q0.front();
                q0.pop();

                if (!currentLevel)
                    continue;

                if (reverse)
                    G.forInNeighborsOf(u, processNeighbor);
                else
                    G.forNeighborsOf(u, processNeighbor);

            } while (!q0.empty());

            std::swap(q0, q1);
        } while (!q0.empty() && currentLevel < bfsLevel);
    };

    kHopGraph.parallelForNodes([&](node u) { prunedBFS(u, K, false); });

    if (kHopGraph.isDirected())
        kHopGraph.parallelForNodes([&](node u) { prunedBFS(u, K, true); });
}

} /* namespace NetworKit */
