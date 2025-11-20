/*
 * EdgeScore.hpp
 *
 *  Created on: 18.08.2015
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_EDGESCORES_EDGE_SCORE_HPP_
#define NETWORKIT_EDGESCORES_EDGE_SCORE_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
/**
 * Abstract base class for an edge score.
 */
template <typename T>
class EdgeScore : public Algorithm {

public:
    EdgeScore(const Graph &G);

    /** Compute the edge score. */
    void run() override;

    /** Get a vector containing the score for each edge in the graph.
     *
     * @return the edge scores calculated by @ref run().
     */
    const std::vector<T> &scores() const;

    /** Get the edge score of the edge with the given edge id.
     */
    T score(edgeid eid);

    /** Get the edge score of the given edge.
     */
    T score(node u, node v);

    /**
     * Build a weighted graph from the computed edge scores.
     *
     * The returned graph has the same topology as the input graph but
     * edge weights are derived from the edge scores:
     *
     *   w(e) = offset + factor * score(e)
     *   or, if squared == true:
     *   w(e) = offset + factor * score(e)^2
     *
     * Requires that edges are indexed (G->hasEdgeIds()).
     * Currently intended for undirected graphs (same behavior as the
     * former EdgeScoreAsWeight helper).
     */
    Graph calculate(bool squared = false, edgeweight offset = 1, edgeweight factor = 1) const;

protected:
    const Graph *G;
    std::vector<T> scoreData;
};

} // namespace NetworKit

#endif // NETWORKIT_EDGESCORES_EDGE_SCORE_HPP_
