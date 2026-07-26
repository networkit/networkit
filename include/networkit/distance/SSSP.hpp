/*
 * SSSP.hpp
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#ifndef NETWORKIT_DISTANCE_SSSP_HPP_
#define NETWORKIT_DISTANCE_SSSP_HPP_

#include <set>

#include <networkit/auxiliary/Multiprecision.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Abstract base class for single-source shortest path algorithms.
 */
template <class GraphT>
class SingleSourceShortestPaths : public Algorithm {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    static constexpr NodeT nullNodeId = NullNodeId<NodeT>;

public:
    /**
     * Creates the SingleSourceShortestPaths class for @a G and source @a s.
     *
     * @param G The graph.
     * @param source The source node.
     * @param storePaths Paths are reconstructable and the number of paths is
     * stored.
     * @param storeNodesSortedByDistance Store a vector of nodes ordered in
     * increasing distance from the source.
     * @param target The target node.
     */
    SingleSourceShortestPaths(const GraphT &G, NodeT source, bool storePaths = true,
                              bool storeNodesSortedByDistance = false, NodeT target = nullNodeId);

    ~SingleSourceShortestPaths() override = default;

    /** Computes the shortest paths from the source to all other nodes. */
    void run() override = 0;

    /**
     * Returns a vector of weighted distances from the source node, i.e. the
     * length of the shortest path from the source node to any other node.
     *
     * @return The weighted distances from the source node to any other node in
     * the graph.
     */
    const std::vector<EdgeWeightT> &getDistances() { return distances; }

    /**
     * Returns the distance from the source node to @a t.
     * @param  t Target node.
     * @return The distance from source to target node @a t.
     */
    EdgeWeightT distance(NodeT t) const { return distances[t]; }

    /**
     * Returns the number of shortest paths between the source node and @a t.
     * @param  t Target node.
     * @return The number of shortest paths between source and @a t.
     */
    bigfloat numberOfPaths(NodeT t) const {
        if (!storePaths) {
            throw std::runtime_error("number of paths have not been stored");
        }
        return npaths[t];
    }

    /**
     * Returns the number of shortest paths between the source node and @a t
     * as a double value. Workaround for Cython
     * @param  t Target node.
     * @return The number of shortest paths between source and @a t.
     */
    double _numberOfPaths(NodeT t) const;

    /**
     * Returns the predecessor nodes of @a t on all shortest paths from source
     * to @a t.
     * @param t Target node.
     * @return The predecessors of @a t on all shortest paths from source to @a
     * t.
     */
    const std::vector<NodeT> &getPredecessors(NodeT t) const {
        if (!storePaths) {
            throw std::runtime_error("predecessors have not been stored");
        }
        return previous[t];
    }

    /**
     * Returns a shortest path from source to @a t and an empty path if source
     * and @a t are not connected.
     *
     * @param t Target node.
     * @param forward If @c true (default) the path is directed from source to
     * @a t, otherwise the path is reversed.
     * @return A shortest path from source to @a t or an empty path.
     */
    std::vector<NodeT> getPath(NodeT t, bool forward = true) const;

    /**
     * Returns all shortest paths from source to @a t and an empty set if source
     * and @a t are not connected.
     *
     * @param t Target node.
     * @param forward If @c true (default) the path is directed from source to
     * @a t, otherwise the path is reversed.
     * @return All shortest paths from source node to target node @a t.
     */
    std::set<std::vector<NodeT>> getPaths(NodeT t, bool forward = true) const;

    /* Returns the number of shortest paths to node t.*/
    bigfloat getNumberOfPaths(NodeT t) const { return npaths[t]; }

    /**
     * Returns a vector of nodes ordered in increasing distance from the source.
     *
     * For this functionality to be available, storeNodesSortedByDistance has
     * to be set to true in the constructor. There are no guarantees regarding
     * the ordering of two nodes with the same distance to the source.
     *
     * @return vector of nodes ordered in increasing distance from the source
     */
    const std::vector<NodeT> &getNodesSortedByDistance() const;

    /**
     * Returns the number of nodes reached by the source.
     *
     * @return Number of nodes reached by the source.
     */
    count getReachableNodes() const {
        assureFinished();
        return reachedNodes;
    }

    /**
     * Sets a new source.
     *
     * @param newSource The new source node.
     */
    void setSource(NodeT newSource) {
        if (!G->hasNode(newSource))
            throw std::runtime_error("Error: node not in the graph.");
        source = newSource;
    }

    /**
     * Sets a new target.
     */
    void setTarget(NodeT newTarget) {
        if (!G->hasNode(newTarget))
            throw std::runtime_error("Error: node not in the graph.");
        target = newTarget;
    }

    /**
     * Returns the sum of distances from the source node to the reached
     * nodes.
     */
    double getSumOfDistances() const {
        assureFinished();
        return sumDist;
    }

protected:
    const GraphT *G;
    NodeT source;
    NodeT target;
    double sumDist;
    count reachedNodes;
    std::vector<EdgeWeightT> distances;
    std::vector<std::vector<NodeT>> previous; // predecessors on shortest path
    std::vector<bigfloat> npaths;

    std::vector<NodeT> nodesSortedByDistance;

    bool storePaths;                 //!< if true, paths are reconstructable and the number of
                                     //!< paths is stored
    bool storeNodesSortedByDistance; //!< if true, store a vector of nodes
                                     //!< ordered in increasing distance from
                                     //!< the source
};

template <class GraphT>
inline double SingleSourceShortestPaths<GraphT>::_numberOfPaths(NodeT t) const {
    if (!storePaths) {
        throw std::runtime_error("number of paths have not been stored");
    }
    bigfloat limit = std::numeric_limits<double>::max();
    if (npaths[t] > limit) {
        throw std::overflow_error("number of paths do not fit into a double");
    }
    double res;
    npaths[t].ToDouble(res);
    return res;
}

using SSSP = SingleSourceShortestPaths<AdjListGraph<node, edgeweight>>;

} /* namespace NetworKit */

#include <networkit/distance/SSSPImpl.hpp>

#endif // NETWORKIT_DISTANCE_SSSP_HPP_
