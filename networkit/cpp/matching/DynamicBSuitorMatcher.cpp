/*
 * DynamicBSuitorMatcher.cpp
 *
 *  Created on: 06.01.2025
 *      Author: Fabian Brandt-Tumescheit
 *              Frieda Gerharz
 */

#include <cassert>
#include <chrono>

#include <networkit/matching/DynamicBSuitorMatcher.hpp>

namespace NetworKit {

void DynamicBSuitorMatcher::processEdgeInsertion(const WeightedEdge &edge) {

    node u = edge.u;
    node v = edge.v;
    edgeweight w = G->weight(u, v);

    MatchingNode startU = suitors[u].insert({v, w});
    MatchingNode startV = suitors[v].insert({u, w});

    if (startU.id != none) {
        suitors[startU.id].remove(u);
    }

    if (startV.id != none) {
        suitors[startV.id].remove(v);
    }

    if (startU.id != none) {
        trackUpdatePath(0, startU.id);
    }

    if (startV.id != none) {
        trackUpdatePath(0, startV.id);
    }
}

void DynamicBSuitorMatcher::trackUpdatePath(size_t batchId, node start, bool recursiveCall) {
    bool done = false;

    node current = start;
    node partner = suitors[current].min.id;
    auto heaviest = suitors[current].min.weight;
    edgeweight prev = std::numeric_limits<edgeweight>::max();

    std::vector<MatchingNode> looseEnds;

    do {
        done = true;

        G->forNeighborsOf(current, [&](node x, edgeweight weight) {
            if (suitors[current].hasPartner(x))
                return;

            const auto z = suitors[x].min;

            if ((weight > heaviest || (weight == heaviest && x < partner))
                && (weight > z.weight || (weight == z.weight && current < z.id))
                && (weight <= prev)) {
                partner = x;
                heaviest = weight;
                return;
            }
        });

        if (partner == none || heaviest < suitors[partner].min.weight
            || (heaviest == suitors[partner].min.weight && current > suitors[partner].min.id)) {
            break;
        }

        MatchingNode prevCurrent = suitors[current].insert({partner, heaviest});
        MatchingNode prevPartner = suitors[partner].insert({current, heaviest});

        if (prevCurrent.id != none) {
            suitors[prevCurrent.id].remove(current);
            looseEnds.emplace_back(
                MatchingNode{prevCurrent.id, suitors[prevCurrent.id].min.weight});
        }

        if (prevPartner.id != none) {
            suitors[prevPartner.id].remove(partner);
            current = prevPartner.id;
            done = false;
        }

        prev = heaviest;
        partner = suitors[current].min.id;
        heaviest = suitors[current].min.weight;
    } while (!done);

    for (auto &looseEnd : looseEnds) {

        trackUpdatePath(batchId, looseEnd.id, true);
    }
}

void DynamicBSuitorMatcher::processEdgeRemoval(const Edge &edge) {

    node u = edge.u;
    node v = edge.v;

    suitors[u].remove(v);
    suitors[v].remove(u);

    trackUpdatePath(0, u);
    trackUpdatePath(0, v);
}

void DynamicBSuitorMatcher::addEdges(std::vector<WeightedEdge> &edges, bool sort) {
    if (sort) {
        std::sort(edges.begin(), edges.end(),
                  [](const WeightedEdge &a, const WeightedEdge &b) { return a.weight > b.weight; });
    }

    for (const auto &edge : edges) {
        if ((suitors[edge.u].hasPartner(edge.v) && suitors[edge.v].hasPartner(edge.u))
            || !isBetterMatch(edge.u, edge.v, edge.weight)
            || !isBetterMatch(edge.v, edge.u, edge.weight)) {
            continue;
        }

        processEdgeInsertion(edge);
    }
}

void DynamicBSuitorMatcher::removeEdges(std::vector<Edge> &edges) {
    for (const auto &edge : edges) {
        assert(!G->hasEdge(edge.u, edge.v));
        if (suitors[edge.u].hasPartner(edge.v)) {

            processEdgeRemoval(edge);
        }
    }
}

} // namespace NetworKit
