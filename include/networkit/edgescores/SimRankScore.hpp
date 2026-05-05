/*  SimRankScore.hpp
 *
 *  Created on: 01.05.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_EDGESCORES_SIM_RANK_SCORE_HPP_
#define NETWORKIT_EDGESCORES_SIM_RANK_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Computes a SimRank similarity score for each edge of the input graph.
 *
 * The score assigned to an edge (u, v) is the SimRank similarity between
 * its endpoints u and v.
 *
 * This implementation follows the classical iterative SimRank recurrence:
 * two nodes are similar if their neighbors/predecessors are similar.
 *
 * The result is stored as one score per edge id in scoreData.
 *
 * Requires indexed edges, i.e. G.indexEdges() must have been called before run().
 */
class SimRankScore final : public EdgeScore<double> {
public:
    /**
     * Constructs a SimRank edge-score algorithm.
     *
     * @param G The input graph.
     * @param similarityPropagationFactor Factor in [0, 1] that scales the contribution
     *        of neighboring node-pair similarities in each SimRank iteration.
     *        Larger values preserve more propagated similarity.
     * @param maxIterations Maximum number of iterations.
     * @param tolerance Convergence tolerance for the maximum absolute score change
     *        between two consecutive iterations.
     *
     * @throws std::invalid_argument if similarityPropagationFactor is not in [0, 1].
     * @throws std::invalid_argument if maxIterations is 0.
     * @throws std::invalid_argument if tolerance is negative.
     */
    SimRankScore(const Graph &G, double similarityPropagationFactor = 0.9,
                 count maxIterations = 100, double tolerance = 1e-4);
    /**
     * Computes one SimRank score per indexed edge.
     *
     * After run(), scores can be accessed through:
     * - scores()
     * - score(edgeid)
     * - score(u, v)
     *
     * @throws std::runtime_error if the graph has no indexed edges.
     *         Call G.indexEdges() before running the algorithm.
     * @throws std::overflow_error if the internal pairwise score matrix would
     *         exceed the supported index range.
     */
    void run() override;

private:
    double similarityPropagationFactor;
    count maxIterations;
    double tolerance;
};

} // namespace NetworKit
#endif // NETWORKIT_EDGESCORES_SIM_RANK_SCORE_HPP_
