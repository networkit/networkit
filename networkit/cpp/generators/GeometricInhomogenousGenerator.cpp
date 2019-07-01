/*
 * GeometricInhomogenousGenerator.cpp
 *
 *  Created on: 09.05.2019
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <networkit/generators/GeometricInhomogenousGenerator.hpp>
#include "girgs/Generator.hpp"

namespace NetworKit {


GeometricInhomogenousGenerator::GeometricInhomogenousGenerator(count n, double avgDegree, double exp, double alpha, unsigned dim) :
    alpha(alpha)
    // do not initialize pointWeights and pointPosition here to check paremeters before
{
    if (alpha <= 1)
        throw std::runtime_error("Temperature T has to satiesfy 0 <= T < 1");

    if (!(1 <= dim && dim <= 5))
        throw std::runtime_error("Support only 1 to 5 dimensions");

    if (n >= static_cast<count>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Support only INT_MAX nodes");

    pointPositions = girgs::generatePositions(static_cast<int>(n), dim, true);
    pointWeights = girgs::generateWeights(static_cast<int>(n), exp, true);
    girgs::scaleWeights(pointWeights, avgDegree, dim, alpha);
}

GeometricInhomogenousGenerator::GeometricInhomogenousGenerator(std::vector<Coordinate> points, std::vector<double> weights, double avgDegree, double alpha) :
    alpha(alpha),
    pointPositions(std::move(points)),
    pointWeights(std::move(weights))
{
    checkInputParameters();
    const auto dim = static_cast<int>(pointPositions.front().size());
    girgs::scaleWeights(pointWeights, avgDegree, dim, alpha);
}

// Construct without scaling
GeometricInhomogenousGenerator::GeometricInhomogenousGenerator(std::vector<Coordinate> points, std::vector<double> weights, double alpha) :
    alpha(alpha),
    pointPositions(std::move(points)),
    pointWeights(std::move(weights))
{
    checkInputParameters();
}

Graph GeometricInhomogenousGenerator::generate() {
    return girgs::generateEdges(pointWeights, pointPositions, alpha, false);
}

Graph GeometricInhomogenousGenerator::generateKeepingInput() {
    return girgs::generateEdges(pointWeights, pointPositions, alpha, true);
}

void GeometricInhomogenousGenerator::checkInputParameters() const {
    if (alpha <= 1)
        throw std::runtime_error("Alpha has to be larger than 1");


    if (!pointPositions.size())
        throw std::runtime_error("Generator requires at least one point");


    if (pointPositions.size() != pointWeights.size())
        throw std::runtime_error("Number of nodes in vectors points and weight do not match");


    const auto dim = pointPositions.front().size();
    if (!(1 <= dim && dim <= 5))
        throw std::runtime_error("Support only 1 to 5 dimensions");


    #ifndef NDEBUG
    for(const auto& pt : pointPositions) {
        if (pt.size() != dim)
            throw std::runtime_error("All points have to have the same dimension");


        for(const auto x : pt) {
            if (!(0 <= x && x <= 1))
                throw std::runtime_error("Points have to lie within the [0:1]^d unit hypercube");

        }
    }
    #endif
}

} // NetworKit
