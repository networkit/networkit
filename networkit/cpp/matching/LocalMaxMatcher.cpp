/*
 * LocalMaxMatcher.cpp
 *
 *  Created on: 05.12.2012
 */

#include <stdexcept>

#include <networkit/matching/LocalMaxMatcher.hpp>

namespace NetworKit {

LocalMaxMatcher::LocalMaxMatcher(const Graph &G) : Matcher(G) {
    if (G.isDirected())
        throw std::runtime_error("Matcher only defined for undirected graphs");
}

// TODO: update to new edge attribute system
// TODO: make local max matching parallel

void LocalMaxMatcher::run() {
    count z = G->upperNodeIdBound();
    count E = G->numberOfEdges();

    std::vector<WeightedEdge> edges;
    edges.reserve(G->numberOfEdges());
    for (auto edge : G->edgeWeightRange())
        edges.emplace_back(edge.u, edge.v, edge.weight + Aux::Random::real(1e-6));

    // candidates records mating candidates
    std::vector<WeightedEdge> candidates(z);
    G->parallelForNodes([&](node u) {
        candidates[u].weight = (edgeweight)0;
        candidates[u].v = u; // itself as mating partner => unmatched
    });

    while (E > 0) {
        // for each edge find out if it is locally maximum
        for (const auto &edge : edges) {
            if (edge.weight > candidates[edge.u].weight
                && edge.weight > candidates[edge.v].weight) {
                candidates[edge.u].v = edge.v;
                candidates[edge.u].weight = edge.weight;
                candidates[edge.v].v = edge.u;
                candidates[edge.v].weight = edge.weight;
            }
        }

        // check if candidates agree to match; if so, then match them
        for (const auto &edge : edges) {
            node u = edge.u;
            node v = edge.v;
            if (candidates[u].v == v && candidates[v].v == u && u != v) {
                // both nodes agree
                M.match(u, v);
            }
        }

        // create remaining "graph" by selecting remaining edges (as triples)
        // adjust candidates
        std::vector<WeightedEdge> newEdges;
        for (const auto &edge : edges) {
            if (!M.isMatched(edge.u) && !M.isMatched(edge.v) && edge.u != edge.v) {
                newEdges.push_back(edge);
                candidates[edge.u].weight = (edgeweight)0;
                candidates[edge.v].weight = (edgeweight)0;
            }
        }
        edges = newEdges;
        E = edges.size();
    }

    hasRun = true;
}

} /* namespace NetworKit */
