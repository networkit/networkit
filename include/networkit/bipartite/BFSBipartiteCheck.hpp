/*  BFSBipartiteCheck.hpp
 *
 *	Created on: 09.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#ifndef NETWORKIT_BFS_BIPARTITE_BIPARTITE_CHECK_HPP_
#define NETWORKIT_BFS_BIPARTITE_BIPARTITE_CHECK_HPP_
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
namespace NetworKit {
class BFSBipartiteCheck final : public Algorithm {

public:
    /**
     * Implements a BFS-based algorithm to check whether a given graph is bipartite.
     * A graph is bipartite if its vertices can be divided into two disjoint sets such that
     * no two adjacent vertices belong to the same set.
     *
     * This algorithm uses a Breadth-First Search (BFS) traversal to attempt a two-coloring
     * of the graph. If a conflict is found (i.e., two adjacent nodes are assigned the same
     * color), the graph is not bipartite.
     *
     * The algorithm runs in O(V + E) time complexity, where:
     * - V is the number of vertices.
     * - E is the number of edges.
     *
     * @param G The input graph to check for bipartiteness. The graph should be undirected.
     * @throws std::runtime_error if the input graph is directed, as bipartiteness is
     *         typically defined for undirected graphs.
     */
    BFSBipartiteCheck(const Graph &G) : graph(&G) {
        if (G.isDirected()) {
            throw std::runtime_error("The graph is not an undirected graph.");
        }
    }

    /**
     * Executes the BFS-based bipartiteness check on the input graph.
     * This method performs all necessary computations to determine whether
     * the graph is bipartite by attempting to two-color the vertices.
     *
     * The result of the bipartiteness check can be retrieved using the `isBipartite()` method.
     */
    void run() override;

    /**
     * Returns whether the input graph is bipartite.
     * The result is only valid after the `run()` method has been called.
     *
     * @return True if the graph is bipartite, false otherwise.
     * @throws std::runtime_error if called before `run()` has been executed.
     */
    bool isBipartite() const {
        assureFinished();
        return isGraphBiPartite;
    }

private:
    const Graph *graph;
    bool isGraphBiPartite{};
};
} // namespace NetworKit
#endif // NETWORKIT_BFS_BIPARTITE_BIPARTITE_CHECK_HPP_
