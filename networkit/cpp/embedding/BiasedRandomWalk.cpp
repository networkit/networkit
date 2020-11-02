/*
 * BiasedRandomWalk.cpp
 *
 *  Created on: 03.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented in C++17
 *  from node2vec [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/Random.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"

namespace NetworKit {
namespace Embedding {

using NeighborSet = std::unordered_set<node>;
using NeighborMap = std::unordered_map<node, AliasSampler>;

struct GraphData {
    // assuming contiguous node numbers 0..n-1, they serve as indices
    using Data = std::vector<NeighborMap>;
    Data data;

    GraphData(size_t nn) : data{nn, NeighborMap()} 
    {}
};

std::unique_ptr<GraphData> graphData;

void preprocessNode(const Graph &graph, node t, double paramP, double paramQ) {

    // for node t
    NeighborSet tNbrs; // Neighbors of t
    for (auto i : Graph::NeighborRange<>(graph, t)) {
        tNbrs.insert(i);
    }
    // for each t-neighbor v:
    for (auto v : Graph::NeighborRange<>(graph, t)) {
        double pSum = 0;
        std::vector<float> pTable; // Probability distribution table
        // for each v-neighbor x:
        for (auto x : Graph::NeighborRange<>(graph, v)) {
            auto weight = graph.weight(v, x);
            if (x == t) {
                pTable.push_back(weight / paramP);
                pSum += weight / paramP;
            } else if (tNbrs.count(x)) {
                pTable.push_back(weight);
                pSum += weight;
            } else {
                pTable.push_back(weight / paramQ);
                pSum += weight / paramQ;
            }
        }
        // Normalizing table
        float pfSum = (float)pSum;
        std::for_each(pTable.begin(), pTable.end(), [pfSum](float &p) { p /= pfSum; });

        graphData->data[v][t].unigram(pTable);
    }
}

std::vector<std::vector<node>> index2node;

// Preprocess transition probabilities for each path t->v->x
void preprocessTransitionProbs(
    const Graph &graph, 
    double paramP, 
    double paramQ 
    ) {

    auto nn = graph.numberOfNodes();
    //graphData = std::make_unique<GraphData>(nn);
    // back to c++11 ;-(
    graphData = std::unique_ptr<GraphData>(new GraphData(nn));
    
    index2node.resize(nn);

    // pre-allocate unordered maps:
    for (index v = 0; v < nn; ++v) {
        auto degv = graph.degreeOut(v);
        auto neighbors = Graph::NeighborRange<>(graph, v);
        for (auto n: neighbors) {
            // init index2node: 
            index2node[v].push_back(n);
            graphData->data[v][n] = AliasSampler(degv);
        }
    }

    count nCnt = 0;

#pragma omp parallel for schedule(dynamic)
    // is forNodes parallelized for omp ?
//  for (node t = 0; t < nn; ++t) {
    for (omp_index t = 0; t < static_cast<omp_index>(nn); ++t) { 
        preprocessNode(graph, t, paramP, paramQ);
        ++nCnt;
    }
}

node nthNeighbor(Graph::NeighborRange<> range,  index nth);


// Simulates a random walk
Walk oneWalk(const Graph &graph, node start, count walkLen) {
    Walk walk(walkLen);
    count nr = 0;
    walk[nr++] = start;
    node src = start;

    if (walkLen == 1) {
        return walk; 
    }
    if (graph.degreeOut(start) == 0) {
        walk.resize(1);       // shorten walk to 1
        return walk;
    }
    auto neighbors = Graph::NeighborRange<>(graph, start);
    auto nn = std::distance(neighbors.begin(), neighbors.end());

    //node randomNeighbor = nthNeighbor(neighbors, nn * uniform_real()); 
    node randomNeighbor = index2node[start][nn * uniform_real()]; 

    walk[nr++] = randomNeighbor;
    node dst   = randomNeighbor;

    while (nr < walkLen) {
        if (graph.degreeOut(dst) == 0) {
            walk.resize(nr);  // shorten walk to nr
            return walk;
        }
        //auto neighbors = Graph::NeighborRange<>(graph, dst);
        //auto nn = std::distance(neighbors.begin(), neighbors.end());
        NeighborMap& map = graphData->data[dst];
        AliasSampler& as = map[src];
        //node next = nthNeighbor(neighbors, as.sample());
        node next = index2node[dst][as.sample()];
        walk[nr++] = next;
        src = dst;
        dst = next;
    }
    return walk;
}

struct WalkData {
    // assuming contiguous node numbers 0..n-1, they serve as indices
    AllWalks data;
    count walkLength;
    count walksPerNode;

    WalkData(count nn, count wl, count wpn) : data{nn * wpn, Walk(wl)}, walkLength(wl), walksPerNode(wpn)
    {}
};

///Simulates walks from every node and writes it into Walks vector
AllWalks doWalks(
    const Graph& graph, 
    count walkLen, 
    count walksPerNode
    ) {
    auto nn = graph.numberOfNodes();
    auto walksDone = 0;
    
    WalkData walkData(nn, walkLen, walksPerNode);

    std::vector<node> shuffled(nn);    
//  std::iota(shuffled.begin(), shuffled.end(), 0);
    int n = 0;
    std::generate(shuffled.begin(), shuffled.end(), [&n](){return n++; });

    for (index c = 0; c < walksPerNode; ++c) {
        std::shuffle(shuffled.begin(), shuffled.end(), getURNG());
#pragma omp parallel for schedule(dynamic)
        // for (index i = 0; i < nn; ++i) { OMP needs signed index var
        for (omp_index i = 0; i < static_cast<omp_index>(nn); ++i) {

            auto v = shuffled[i];
            Walk thisWalk = oneWalk(graph, v, walkLen);            
            walkData.data[c * nn + i] = std::move(thisWalk);

            ++walksDone;
        }
    }
    
    graphData.reset();

    return walkData.data;
}

node nthNeighbor(Graph::NeighborRange<> range,  index nth) {
    auto it = range.begin();
    std::advance(it, nth);
    return *it;
}

} // namespace Embedding
} // namespace NetworKit
