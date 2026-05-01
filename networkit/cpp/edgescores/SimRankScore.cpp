/*  SimRankScore.cpp
*
 *  Created on: 01.05.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/edgescores/SimRankScore.hpp>

namespace NetworKit {

SimRankScore::SimRankScore(const Graph &G, double similarityPropagationFactor, count maxIterations, double tolerance)
    : EdgeScore<double>(G), similarityPropagationFactor(similarityPropagationFactor), maxIterations(maxIterations), tolerance(tolerance), iterations(0)
{
    if (similarityPropagationFactor < 0.0 || similarityPropagationFactor > 1.0) {
        throw std::invalid_argument("damping must be in the range [0,1]");
    }

    if (maxIterations == 0) {
        throw std::invalid_argument("maxIterations must be greater than 0");
    }

    if (tolerance < 0.0) {
        throw std::invalid_argument("tolerance must be greater than or equal to 0");
    }
}

void SimRankScore::run() {
    hasFinished();
}

} // namespace NetworKit
