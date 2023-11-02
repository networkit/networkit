/*
 * GeometricInhomogeneousGenerator.cpp
 *
 *  Created on: 26.10.2023
 *      Author: Christopher Weyand, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <limits>

#include <girgs/Generator.h>

#include <networkit/generators/GeometricInhomogeneousGenerator.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

GeometricInhomogeneousGenerator::GeometricInhomogeneousGenerator(count n, double avgDegree, double exp, double temperature, unsigned dim)
    : n(n), avgDegree(avgDegree), powerlawExp(exp), dim(dim) {
    if (n >= static_cast<count>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Support only INT_MAX nodes");

    if (avgDegree > n-1)
        throw std::runtime_error("Average degree is at most n-1");

    if (exp <= 2.0)
        throw std::runtime_error("Power law exponent must be > 2");

    if (temperature < 0 || temperature > 1)
        throw std::runtime_error("Temperature T has to satiesfy 0 <= T <= 1");
    alpha = temperature <= 0.0 ? std::numeric_limits<double>::infinity() : 1/temperature;

    if (!(1 <= dim && dim <= 5))
        throw std::runtime_error("Support only 1 to 5 dimensions");
}

Graph GeometricInhomogeneousGenerator::generate() {
    auto& gen = Aux::Random::getURNG();

    // generate inputs
    auto positions = girgs::generatePositions(static_cast<int>(n), dim, gen());
    auto weights = girgs::generateWeights(static_cast<int>(n), powerlawExp, gen());
    girgs::scaleWeights(weights, avgDegree, dim, alpha);

    // gerenate edge list
    auto edges = girgs::generateEdges(weights, positions, alpha, gen());

    // convert to networkit Graph
    Graph G(n);
    for(auto [u,v] : edges) 
        G.addEdge(u,v);

    return G;
}

} // NetworKit