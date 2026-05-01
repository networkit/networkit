/*  SimRankScore.cpp
*
 *  Created on: 01.05.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/edgescores/SimRankScore.hpp>

namespace NetworKit {

SimRankScore::SimRankScore(const Graph &G, double damping, count maxIterations, double tolerance)
    : EdgeScore<double>(G), damping(damping), maxIterations(maxIterations), tolerance(tolerance), iterations(0)
{
    if (damping < 0.0 || damping > 1.0) {
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
