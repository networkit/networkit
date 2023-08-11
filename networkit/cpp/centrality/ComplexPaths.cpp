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
#include <networkit/graph/Graph.hpp>
//#include <networkit/graph/GraphTools.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/centrality/ComplexPaths.hpp>

using std::vector;
using std::sample;

namespace NetworKit {

ComplexPathAlgorithm::ComplexPathAlgorithm(const Graph& G, count threshold, Mode mode, node start)
: inputGraph{G}, mode{mode}, 
  start{start}, threshold{threshold}, normPaths{false} {
    if (mode == Mode::singleNode && start == none)
        ERROR("complexPathAlgorithm: provide node != none in Mode::singleNode.");
    if (mode == Mode::allNodes && start != none)
        WARN("complexPathAlgorithm: start node ignored in Mode::allNodes.");
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

vector<double> ComplexPathAlgorithm::getPLci() {
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

vector<node> ComplexPathAlgorithm::getAdopters() {
    if (mode != Mode::singleNode)
        WARN("complexPathAlgorithm[getAdopters]: no results in Mode::allNodes.");
    assureFinished();
    return adopters;
}

static void addNewEdge(Graph& g, node u, node v) {
  if(g.hasEdge(u, v)) return;
  g.addEdge(u,v);
}

Graph ComplexPathAlgorithm::complexPathsGraph(node seed, count threshold, vector<node>* adopters) {
    const auto& g = inputGraph;
    auto n = g.numberOfNodes();
    Graph complex_g = Graph(n); // all nodes, no edges so far
    vector<bool>  activated(n); 
    vector<count> influence(n);

    // activate threshold-1 random neighbors of seed:
    auto needSeeds = threshold - 1;
  
    vector<node> seeds;
    seeds.reserve(4);

    auto neighbors = Graph::NeighborRange<>(g, seed);

    auto randomGen = Aux::Random::getURNG();
    sample(neighbors.begin(), neighbors.end(), 
           std::back_inserter(seeds), needSeeds,  
           randomGen);
    auto needMore = needSeeds - seeds.size();
    if (needMore > 0) {
        WARN("complexPathsGraph: ", seed, ": not enough direct neighbors.");
        vector<node> moreSeeds;
        for (auto u: neighbors) {
            auto indirectNeighbors = Graph::NeighborRange<>(g, u);
            for (auto in: indirectNeighbors) {
            if (in != seed && std::find(moreSeeds.begin(), moreSeeds.end(), in) == moreSeeds.end())
                moreSeeds.push_back(in);
            } 
            sample(moreSeeds.begin(), moreSeeds.end(),
                   std::back_inserter(seeds), needMore, randomGen);  
            if (seeds.size() < needSeeds) 
                WARN("complexPathsGraph: ", seed, ": not enough indirect neighbors.");
        }
    }

    seeds.insert(seeds.begin(), seed); // put in front, so that direct neighbors are activated first
  
    for (auto u: seeds) {
        activated[u] = true;
        if (adopters) adopters->push_back(u);
            auto neighbors = Graph::NeighborRange<>(g, u);

        for (auto v: neighbors) {
            addNewEdge(complex_g, u, v);
        }
    }
    // compute influences: degrees of all nodes in complex_g so far
    auto computeInfluence = [&](node u){ influence[u] = complex_g.degree(u); };
    complex_g.parallelForNodes(computeInfluence);  

    bool spread; 
    do {
        spread = false; // at start of new iteration
        complex_g.forNodes([&](node u) {
            if (influence[u] >= threshold) {
                if (!activated[u]) {
                    activated[u] = true;
                    if (adopters) adopters->push_back(u);
                    spread = true;
                }
                auto uNeighbors = Graph::NeighborRange<>(g, u);
                for (auto v: uNeighbors) {
                    addNewEdge(complex_g, u, v); 
                } 
            } 
        });
        if (spread) { // new influences in growing complex_g
            complex_g.parallelForNodes(computeInfluence);  
        }
    } while(spread);

    return complex_g;
}

static auto min_max_norm(vector<double> data) {
    auto [min, max] = std::minmax_element(data.begin(), data.end());
    auto diff = *max - *min;
    vector<double> norm(data);
    auto scale = [&](auto e) { return (e - *min)/diff; };
    for (auto& e: norm) e = scale(e);

    return norm;
}

vector<double> ComplexPathAlgorithm::complexPathLength(count t) {
    const auto& g = inputGraph;  
    auto n = g.numberOfNodes();
    vector<double> PLci(n);
    constexpr auto infDist = std::numeric_limits<edgeweight>::max();

    g.parallelForNodes([&](node u) {
        auto c_g = complexPathsGraph(u, t, nullptr);
        auto noPaths = false;
        NetworKit::BFS bfs(c_g, u, noPaths);
        bfs.run();
        auto distances = bfs.getDistances();
        for (auto& d: distances) if (d==infDist) d=0.0;
        auto distanceSum = std::accumulate(distances.begin(), distances.end(), 0.0);
        PLci[u] = distanceSum / n;
    });

    if (normPaths) 
    PLci=min_max_norm(PLci);	

    return PLci;
}

} /* namespace NetworKit */

