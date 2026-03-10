/*
 * EpidemicSimulationSEIR.cpp
 *
 *  Created on: 20.11.2015
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/simulation/EpidemicSimulationSEIR.hpp>

namespace NetworKit {

EpidemicSimulationSEIR::EpidemicSimulationSEIR(const Graph &G, count tMax, double transP,
                                               count eTime, count iTime, node start)
    : Algorithm(), G(&G), maxTimestamps(tMax), transmissionProbability(transP), exposedTime(eTime),
      infectiousTime(iTime), start(start) {}

void EpidemicSimulationSEIR::run() {

    using State = EpidemicSimulationSEIR::State;

    index t = 0;

    // initialize state and timestamp arrays
    state.resize(G->upperNodeIdBound(), State::Undefined);
    timestamp.resize(G->upperNodeIdBound(), none);

    auto setState = [&](node v, State X) {
        state[v] = X;
        timestamp[v] = t;
    };

    // initialize nodes to Susceptible
    G->parallelForNodes([&](node v) { setState(v, State::Susceptible); });

    // contact may expose susceptible node to infection
    auto contact = [&](node v) {
        if ((state[v] == State::Susceptible)
            && (Aux::Random::probability() <= transmissionProbability)) {
            setState(v, State::Exposed);
        }
    };

    // update state of nodes
    auto sweep = [&](node u) {
        if (state[u] == State::Susceptible) {
            // do nothing
        } else if (state[u] == State::Exposed) {
            // exposed nodes become infectious after time
            if ((t - timestamp[u]) >= exposedTime) {
                setState(u, State::Infectious);
            }
        } else if (state[u] == State::Infectious) {
            // contact neighbors of infectious node
            G->forNeighborsOf(u, [&](node v) { contact(v); });
            // infectious nodes become removed after time
            if ((t - timestamp[u]) >= infectiousTime) {
                setState(u, State::Removed);
            }
        } else if (state[u] == State::Removed) {
            // do nothing
        } else if (state[u] == State::Undefined) {
            throw std::runtime_error("node in undefined state encountered - should not happen");
        } else {
            throw std::runtime_error("else branch taken - should not happen");
        }
    };

    auto census = [&]() {
        std::vector<count> data(5);
        G->forNodes([&](node v) { data[(index)state[v]] += 1; });
        return data;
    };

    // if starting node node provided, start with random node
    if (start == none) {
        start = GraphTools::randomNode(*G);
    }
    INFO("zero node: ", start);
    setState(start, State::Infectious); // infect node zero

    while (t < maxTimestamps) {
        G->parallelForNodes(sweep);
        auto populations = census();

        for (int s = (int)State::Susceptible; s != (int)State::Undefined; ++s) {
            std::vector<count> data = {start, t, (count)s, populations[s]};
            stats.push_back(data);
        }

        t += 1;
    }

    hasRun = true;
}

const std::vector<std::vector<count>> &EpidemicSimulationSEIR::getData() const {
    return stats;
}

} /* namespace NetworKit */
