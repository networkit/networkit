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
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"
#include "LearnEmbedding.hpp"
#include <networkit/embedding/KHop.hpp>

namespace NetworKit {

KHop::KHop(const Graph &G, size_t K, double S, count L, count N, count D, khopMode M, count winSize,
           count iterations)
    : G(&G), K(K), S(S), L(L), N(N), D(D), M(M), winSize(winSize), iterations(iterations) {

    if (S < 0.0)
        throw std::runtime_error("S should be positive.");
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

    // Generate G_k

    // create G with same set of nodes
    this->G_k = NetworKit::Graph(this->G->numberOfNodes(), true, this->G->isDirected());
    // construct G_k by drawing edges to k-hop neighbors
    this->G->forNodes([&](node v) { this->kHop(v); });

    // set probabilities and assign them as edge-weights
    double min = this->G_k.numberOfNodes();
    double max = 0.0;
    this->G_k.forNodes([&](node v) {
        auto begin = this->G_k.neighborRange(v).begin();
        auto end = this->G_k.neighborRange(v).end();

        // if the neighbors of v have a distance of 0, self-loop is the only edge of v
        if (std::distance(begin, end) == 0) {
            this->G_k.addEdge(v, v); // edged on weighted graphs have a default weight of 1.0
        } else {
            this->G_k.forNeighborsOf(v, [&](node u) {
                double weight = std::exp(this->overlapCoefficient(v, u));

                // assign computed weight
                this->G_k.setWeight(v, u, weight);

                if (weight < min)
                    min = weight;
                if (weight > max)
                    max = weight;
            });
        }
    });

    // scale weights
    this->G_k.forEdges([&](node v, node u, edgeweight w) {
        this->G_k.setWeight(v, u, this->rescale(w, min, max, 1.0, 1 + this->S));
    });

    // normalize
    this->G_k.forNodes([&](node v) {
        double sum = 0;
        count degree = this->G_k.degreeOut(v);
        std::vector<double> weightsByNeighborIndex(degree, 0.0);

        // get number of neighbors and iterate (faster than edgeweight weight(node u, node v)
        // const;)
        for (count i = 0; i < degree; ++i) {
            weightsByNeighborIndex[i] = this->G_k.getIthNeighborWeight(v, i);
            sum += weightsByNeighborIndex[i];
        }

        for (count i = 0; i < degree; ++i) {
            this->G_k.setWeight(v, this->G_k.getIthNeighbor(v, i), weightsByNeighborIndex[i] / sum);
        }
    });
}

void KHop::run() {
    using namespace Embedding;

    Aux::SignalHandler handler;

    TRACE("preprocess transition probabilities ...");

    BiasedRandomWalk brw(&this->G_k);
    brw.preprocessTransitionProbs(1, 1);
    handler.assureRunning();

    TRACE("do biased walks ...");
    auto walks = brw.doWalks(this->L, this->N);
    handler.assureRunning();

    TRACE("learn embeddings ...");
    features =
        learnEmbeddings(walks, this->G_k.numberOfNodes(), this->D, this->winSize, this->iterations);
    handler.assureRunning();

    hasRun = true;
}

const std::vector<std::vector<float>> &KHop::getFeatures() const {
    assureFinished();
    return features;
}

double KHop::overlapCoefficient(node node1, node node2) {
    return (double)NeighborhoodUtility::getCommonNeighbors(*this->G, node1, node2).size()
           / (double)(std::min(this->G->degree(node1), this->G->degree(node2)));
}

double KHop::rescale(double weight, double min, double max, double low, double high) {
    return ((weight - min) / (max - min)) * (high - low) + low;
}

void KHop::kHop(node node) {

    // queue of tuples with nodes and their 'depth'
    std::queue<std::tuple<NetworKit::node, count>> queue;
    // vector of visited nodes and the depth at which they have been visited (while -1 means not
    // visited)
    std::vector<int64_t> visited(this->G->numberOfNodes(), -1);

    visited[node] = true;
    queue.emplace(std::make_tuple(node, 0));

    while (!queue.empty()) {
        std::tuple<NetworKit::node, count> q = queue.front();
        queue.pop();

        NetworKit::node current = std::get<0>(q);
        count depth = std::get<1>(q);

        if (depth < this->K) {
            this->G->forNeighborsOf(current, [&](NetworKit::node neighbor) {
                if (visited[neighbor] == -1 || static_cast<count>(visited[neighbor]) > depth + 1) {
                    visited[neighbor] = depth + 1;
                    queue.push(std::make_tuple(neighbor, depth + 1));
                }
            });
        }

        if (((this->M == khopMode::STRICT && depth == this->K)
             || (this->M == khopMode::DEFAULT && depth <= this->K))
            && !this->G_k.hasEdge(node, current)) {
            this->G_k.addEdge(node, current);
        }
    }
}

} /* namespace NetworKit */
