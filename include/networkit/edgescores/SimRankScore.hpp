/*  SimRankScore.hpp
 *
 *  Created on: 01.05.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_EDGESCORES_SIMRANK_SCORE_HPP
#define NETWORKIT_EDGESCORES_SIMRANK_SCORE_HPP

#include <networkit/edgescores/EdgeScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class SimRankScore final : public EdgeScore<double> {
public:
    SimRankScore(const Graph &G, double damping = 0.9, count maxIterations = 100,
                 double tolerance = 1e-4);

    void run() override;

private:
    double damping;
    count maxIterations;
    double tolerance;

    count iterations;
};

} // namespace NetworKit
#endif // NETWORKIT_EDGESCORES_SIMRANK_SCORE_HPP
