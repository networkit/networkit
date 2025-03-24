/*
 * DYNBETNODE.hpp
 *
 *  Created on: 10.03.2016
 *      Author: Elisabetta Bergamini
 */

#ifndef NETWORKIT_CENTRALITY_DYN_BETWEENNESS_ONE_NODE_HPP_
#define NETWORKIT_CENTRALITY_DYN_BETWEENNESS_ONE_NODE_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * This algorithm computes the betweenness centrality for a specified focus node x.
 * It works for dynamic graphs that can be undirected/directed and unweighted/weighted without
 * negative cycles. The Algorithm proceeds as described in [1]. The run() method computes the
 * betweenness centrality by computing SSSP (similar to Dijkstra) for each node.
 * The update() method only works with edge insertions and edge weight decrements and has a WC
 * runtime complexity of O(n²) which outperforms the DynamicBetweenness algorithm for all nodes
 * (instead of one focus node) that takes O(nm).
 *
 * [1] Bergamini et al. : Improving the Betweenness Centrality of a Node by Adding Links
 *    https://dl.acm.org/doi/abs/10.1145/3166071
 */

class DynBetweennessOneNode : public Algorithm, public DynAlgorithm {

public:
    /**
     * Creates the object for @a G.
     *
     * @param G The graph.
     * @paarm x The node for which we want to compute betweenness.
     */
    DynBetweennessOneNode(Graph &G, node x);

    /** initialize distances and Pred by repeatedly running the Dijkstra2 algorithm */
    void run() override;

    /* Computes the betweenness centrality of x by running BFS for each node */
    void runUnweighted();

    /* Computes the betweenness centrality of x by running Dijkstra for each node */
    void runWeighted();

    /**
     * Updates the betweenness centrality of x after an edge insertions on the graph.
     * Notice: it works only with edge insertions.
     *
     * @param e The edge insertions.
     */
    void update(GraphEvent event) override;

    /**
     * Updates the betweenness centrality of x after a batch of edge insertions on the graph.
     * Notice: it works only with edge insertions.
     *
     * @param batch The batch of edge insertions.
     */
    void updateBatch(const std::vector<GraphEvent> &batch) override;

    /** Returns the new betweenness score of node x after the insertion of an edge. Distances and
     * scores of the other nodes are not changed by this function. */
    edgeweight computeScore(GraphEvent event);

    edgeweight getDistance(node u, node v);
    edgeweight getSigma(node u, node v);
    edgeweight getSigmax(node u, node v);
    edgeweight getbcx();

private:
    Graph &G;
    node x;
    // betweenness centrality of node x
    edgeweight bcx = 0;
    const edgeweight infDist = std::numeric_limits<edgeweight>::max();
    const edgeweight epsilon =
        0.0000000001; // make sure that no legitimate edge weight is below that.

    std::vector<std::vector<edgeweight>> distances;
    std::vector<std::vector<edgeweight>> distancesOld;
    // total number of shortest paths between two nodes
    std::vector<std::vector<edgeweight>> sigma;
    // number of shortest paths between two nodes that go through node x
    std::vector<std::vector<edgeweight>> sigmax;
    std::vector<std::vector<edgeweight>> sigmaxOld;

    std::vector<node> Pred;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_DYN_BETWEENNESS_ONE_NODE_HPP_
