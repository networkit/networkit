/*
 *  ComplexPaths.cpp
 *
 *
 *  Created on:16.06.2023
 *      Author: Klaus Ahrens
 *              <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented from
 *
 *            https://github.com/drguilbe/complexpaths.git
 *
 *  see [ Guilbeault, D., Centola, D. Topological measures for
 *        identifying and predicting the spread of complex contagions.
 *        Nat Commun 12, 4430 (2021).
 *        https://doi.org/10.1038/s41467-021-24704-6 ]
 *
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/centrality/ComplexPaths.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

ComplexPathAlgorithm::ComplexPathAlgorithm(const Graph &G, count threshold, Mode mode, node start)
    : inputGraph{&G}, mode{mode}, start{start}, threshold{threshold}, normPaths{false} {
    if (threshold < 1)
        throw std::runtime_error("complexPathAlgorithm: threshold must be greater than 0.");
    if (mode == Mode::singleNode && start == none)
        throw std::runtime_error("complexPathAlgorithm: provide node != none in Mode::singleNode.");
    if (mode == Mode::allNodes && start != none)
        throw std::runtime_error("complexPathAlgorithm: start node ignored in Mode::allNodes.");
}

void ComplexPathAlgorithm::normalize() {
    if (mode != Mode::allNodes) {
        WARN("complexPathAlgorithm: no normalization in Mode::singleNode.");
        return;
    }
    normPaths = true;
}

void ComplexPathAlgorithm::run() {
    if (mode == Mode::singleNode) {
        complexPathGraph = complexPathsGraph(start, threshold, &adopters);
    } else {
        complexPathsLengths = complexPathLength(threshold);
    }
    hasRun = true;
}

std::vector<double> ComplexPathAlgorithm::getPLci() {
    if (mode != Mode::allNodes)
        WARN("complexPathAlgorithm[getPLci]: no results in Mode::singleNode.");
    assureFinished();
    return complexPathsLengths;
}

Graph ComplexPathAlgorithm::getComplexGraph() {
    if (mode != Mode::singleNode)
        WARN("complexPathAlgorithm[getComplexGraph]: no results in Mode::allNodes.");
    assureFinished();
    return complexPathGraph;
}

std::vector<node> ComplexPathAlgorithm::getAdopters() {
    if (mode != Mode::singleNode)
        WARN("complexPathAlgorithm[getAdopters]: no results in Mode::allNodes.");
    assureFinished();
    return adopters;
}

static void addNewEdge(Graph &g, node u, node v) {
    if (g.hasEdge(u, v))
        return;
    g.addEdge(u, v);
}

std::vector<node> ComplexPathAlgorithm::generateSeeds(node seed, const Graph &g, count threshold) {
    // activate threshold-1 random neighbors of seed:
    const count seedsNeeded = threshold - 1;

    std::vector<node> seeds;
    seeds.reserve(4);

    auto neighbors = Graph::NeighborRange<>(g, seed);

    auto randomGen = Aux::Random::getURNG();
    std::sample(neighbors.begin(), neighbors.end(), std::back_inserter(seeds), seedsNeeded,
                randomGen);

    if (seeds.size() >= seedsNeeded)
        return seeds;
    const count needMore = seedsNeeded - seeds.size();

    if (needMore > 0) {
        WARN("complexPathsGraph: ", seed, ": not enough direct neighbors.");

        // collect all nodes at distance 2 before sampling to establish equal distribution
        std::vector<node> indirNeighbors;

        for (node u : neighbors) {
            for (node in : g.neighborRange(u)) {
                if (in != seed) {
                    indirNeighbors.push_back(in);
                }
            }
        }
        std::sample(indirNeighbors.begin(), indirNeighbors.end(), std::back_inserter(seeds),
                    needMore, randomGen);
    }
    return seeds;
}

Graph ComplexPathAlgorithm::complexPathsGraph(node seed, count threshold,
                                              std::vector<node> *adopters) {
    const Graph *g = inputGraph;
    const auto n = g->numberOfNodes();
    Graph complex_g = Graph(n); // all nodes, no edges so far
    std::vector<bool> activated(n);
    std::vector<count> influence(n);

    std::vector<node> seeds = generateSeeds(seed, *g, threshold);
    seeds.insert(seeds.begin(), seed); // put in front, so that direct neighbors are activated first

    for (node u : seeds) {
        activated[u] = true;
        if (adopters)
            adopters->push_back(u);

        for (node v : g->neighborRange(u)) {
            addNewEdge(complex_g, u, v);
        }
    }
    // compute influences: degrees of all nodes in complex_g so far
    auto computeInfluence = [&](node u) { influence[u] = complex_g.degree(u); };
    complex_g.parallelForNodes(computeInfluence);

    bool spread;
    do {
        spread = false; // at start of new iteration
        complex_g.forNodes([&](node u) {
            if (influence[u] >= threshold) {
                if (!activated[u]) {
                    activated[u] = true;
                    if (adopters)
                        adopters->push_back(u);
                    spread = true;
                }
                for (node v : g->neighborRange(u)) {
                    addNewEdge(complex_g, u, v);
                }
            }
        });
        if (spread) { // new influences in growing complex_g
            complex_g.parallelForNodes(computeInfluence);
        }
    } while (spread);

    return complex_g;
}

static std::vector<double> min_max_norm(std::vector<double> &data) {
    auto [min, max] = std::minmax_element(data.begin(), data.end());
    double diff = *max - *min;
    std::vector<double> norm(data);
    // Hint: Use structured bind "min" as init-capture in order to avoid compiler errors
    // for C++17 and Clang (at least up until C++20).
    std::transform(norm.cbegin(), norm.cend(), norm.begin(),
                   [diff, &min = min](double e) { return (e - *min) / diff; });

    return norm;
}

std::vector<double> ComplexPathAlgorithm::complexPathLength(count t) {
    const Graph *g = inputGraph;
    const auto n = g->numberOfNodes();
    std::vector<double> PLci(n);
    const auto infDist = std::numeric_limits<edgeweight>::max();

    g->parallelForNodes([&](node u) {
        Graph c_g = complexPathsGraph(u, t, nullptr);
        bool noPaths = false;
        BFS bfs(c_g, u, noPaths);
        bfs.run();
        auto distances = bfs.getDistances();
        for (auto &d : distances)
            if (d == infDist)
                d = 0.0;
        auto distanceSum = std::accumulate(distances.begin(), distances.end(), 0.0);
        PLci[u] = distanceSum / static_cast<double>(n);
    });

    if (normPaths)
        PLci = min_max_norm(PLci);

    return PLci;
}

} /* namespace NetworKit */
