/*
 * DynamicBSuitorMatcher.cpp
 *
 *  Created on: 06.01.2025
 *      Author: Fabian Brandt-Tumescheit
 *              Frieda Gerharz
 */

#include <cassert>
#include <vector>

#include <networkit/matching/DynamicBSuitorMatcher.hpp>

namespace NetworKit {

void DynamicBSuitorMatcher::update(GraphEvent e) {
    switch (e.type) {
    case GraphEvent::EDGE_ADDITION:
        addEdge(e);
        break;
    case GraphEvent::EDGE_REMOVAL:
        removeEdge(e);
        break;
    default:
        throw std::runtime_error("Event type not allowed. Edge insertions and removals only.");
    }
}

void DynamicBSuitorMatcher::updateBatch(const std::vector<GraphEvent> &batch) {
    for (GraphEvent e : batch) {
        update(e);
    };
}

void DynamicBSuitorMatcher::processEdgeInsertion(const GraphEvent &event) {

    node u = event.u;
    node v = event.v;
    edgeweight w = event.w;

    MatchingNode startU = suitors[u].insert({v, w});
    MatchingNode startV = suitors[v].insert({u, w});

    if (startU.id != none) {
        suitors[startU.id].remove(u);
    }

    if (startV.id != none) {
        suitors[startV.id].remove(v);
    }

    if (startU.id != none) {
        trackUpdatePath(startU.id);
    }

    if (startV.id != none) {
        trackUpdatePath(startV.id);
    }
}

void DynamicBSuitorMatcher::trackUpdatePath(node start) {
    bool done = false;

    node current = start;
    node partner = suitors[current].min.id;
    edgeweight heaviest = suitors[current].min.weight;
    edgeweight prev = std::numeric_limits<edgeweight>::max();

    std::vector<MatchingNode> looseEnds;

    do {
        done = true;

        G->forNeighborsOf(current, [&](node x, edgeweight weight) {
            if (suitors[current].hasPartner(x))
                return;

            const MatchingNode z = suitors[x].min;

            if ((weight > heaviest || (weight == heaviest && x < partner))
                && (weight > z.weight || (weight == z.weight && current < z.id))
                && (weight <= prev)) {
                partner = x;
                heaviest = weight;
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
            looseEnds.emplace_back(prevCurrent.id, suitors[prevCurrent.id].min.weight);
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
        trackUpdatePath(looseEnd.id);
    }
}

void DynamicBSuitorMatcher::processEdgeRemoval(const GraphEvent &event) {

    node u = event.u;
    node v = event.v;

    suitors[u].remove(v);
    suitors[v].remove(u);

    trackUpdatePath(u);
    trackUpdatePath(v);
}

void DynamicBSuitorMatcher::addEdge(const GraphEvent &event) {
    if ((suitors[event.u].hasPartner(event.v) && suitors[event.v].hasPartner(event.u))
        || !isBetterMatch(event.u, event.v, event.w) || !isBetterMatch(event.v, event.u, event.w)) {
        return;
    }
    processEdgeInsertion(event);
}

void DynamicBSuitorMatcher::removeEdge(const GraphEvent &event) {
    if (suitors[event.u].hasPartner(event.v)) {
        processEdgeRemoval(event);
    }
}

} // namespace NetworKit
