/*
 * GeometricInhomogenousGenerator.cpp
 *
 *  Created on: 09.05.2019
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#include "GeometricInhomogenousGenerator.h"
#include "girgs/Generator.h"

namespace NetworKit {


GeometricInhomogenousGenerator::GeometricInhomogenousGenerator(count n, double avgDegree, double exp, double alpha, unsigned dim) :
    alpha(alpha)
    // do not initialize pointWeights and pointPosition here to check paremeters before
{
    if (alpha <= 1)
        throw std::runtime_error("Temperature T has to satiesfy 0 <= T < 1");

    if (!(1 <= dim && dim <= 5))
        throw std::runtime_error("Support only 1 to 5 dimensions");

    pointPositions = girgs::generatePositions(n, dim, true);
    pointWeights = girgs::generateWeights(n, exp, true);
    girgs::scaleWeights(pointWeights, avgDegree, dim, alpha);
}

GeometricInhomogenousGenerator::GeometricInhomogenousGenerator(std::vector<coordinate_t>& points, std::vector<double>& weights, double avgDegree, double alpha) :
    alpha(alpha),
    pointPositions(std::move(points)),
    pointWeights(std::move(weights))
{
    check_input_parameters_();
    const auto dim = pointPositions.front().size();
    girgs::scaleWeights(pointWeights, avgDegree, dim, alpha);
}

// Construct without scaling
GeometricInhomogenousGenerator::GeometricInhomogenousGenerator(std::vector<coordinate_t>& points, std::vector<double>& weights, double alpha) :
    alpha(alpha),
    pointPositions(std::move(points)),
    pointWeights(std::move(weights))
{
    check_input_parameters_();
}

Graph GeometricInhomogenousGenerator::generate() {
    return girgs::generateEdges(pointWeights, pointPositions, alpha, false);
}

Graph GeometricInhomogenousGenerator::generateKeepingInput() {
    return girgs::generateEdges(pointWeights, pointPositions, alpha, true);
}

void GeometricInhomogenousGenerator::check_input_parameters_() const {
    if (alpha <= 1)
        throw std::runtime_error("Alpha has to be larger than 1");


    if (!pointPositions.size())
        throw std::runtime_error("Generator requires at least one point");


    if (pointPositions.size() != pointWeights.size())
        throw std::runtime_error("Number of nodes in vectors points and weight do not match");


    const auto dim = pointPositions.front().size();
    if (!(1 <= dim && dim <= 5))
        throw std::runtime_error("Support only 1 to 5 dimensions");


    for(const auto& pt : pointPositions) {
        if (pt.size() != dim)
            throw std::runtime_error("All points have to have the same dimension");


        for(double x : pt) {
            if (!(0 <= x && x <= 1))
                throw std::runtime_error("Points have to lie within the [0:1]^d unit hypercube");

        }
    }
}

} // NetworKit
