/*
 * EpidemicSimulationSEIR.hpp
 *
 *  Created on: 20.11.2015
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_SIMULATION_EPIDEMIC_SIMULATION_SEIR_HPP_
#define NETWORKIT_SIMULATION_EPIDEMIC_SIMULATION_SEIR_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup simulation
 * Node centrality index which ranks nodes by their degree.
 * Optional normalization by maximum degree.
 */
class EpidemicSimulationSEIR final : public Algorithm {
public:
    /**
     * Simulates an epidemic spread using the Susceptible-Exposed-Infectious-Removed (SEIR) model.
     * @param G The network.
     * @param tMax Max. number of timesteps.
     * @param transP Transmission probability.
     * @param eTime Exposed time.
     * @param iTime Infectious time.
     * @param zero Starting node.
     */
    EpidemicSimulationSEIR(const Graph &G, count tMax, double transP, count eTime, count iTime,
                           node zero);

    void run() override;

    /**
     * Returns the data from the simulation (only valid after run() was called).
     * @return Vector of vectors, containing "zero", "time", "state" and "count" values for each
     * node.
     */
    const std::vector<std::vector<count>> &getData() const;

private:
    const Graph *G;
    count tMax;
    double transP;
    count eTime;
    count iTime;
    node zero;
    enum class State { S, E, I, R, U }; // Susceptible, Exposed, Infectious, Removed, Undefined
    std::vector<State> state;
    std::vector<index> timestamp;
    std::vector<std::vector<count>> stats;
};

} /* namespace NetworKit */

#endif // NETWORKIT_SIMULATION_EPIDEMIC_SIMULATION_SEIR_HPP_
