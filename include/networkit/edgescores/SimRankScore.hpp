//
// Created by andreas on 01.05.26.
//

#ifndef NETWORKIT_SIMRANKSCORE_H
#define NETWORKIT_SIMRANKSCORE_H

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
#endif // NETWORKIT_SIMRANKSCORE_H
