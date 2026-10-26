/*
 * TopologicalSort.hpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */
#ifndef NETWORKIT_GRAPH_TOPOLOGICAL_SORT_HPP_
#define NETWORKIT_GRAPH_TOPOLOGICAL_SORT_HPP_

#include <optional>
#include <unordered_map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Given a directed graph G, the topology sort algorithm creates one valid topology order of nodes.
 * Undirected graphs are not accepted as input, since a topology sort is a linear ordering of
 * vertices such that for every edge u -> v, node u comes before v in the ordering.
 */
template <typename GraphT>
class GenericTopologicalSort final : public Algorithm {
public:
    using NodeT = typename GraphT::NodeT;
    using NodeIdMapping = std::unordered_map<NodeT, index>;

    /**
     * Initialize the topological sort algorithm by passing an input graph.
     *
     * @param G The input graph.
     */
    GenericTopologicalSort(const GraphT &G);

    /**
     * Initialize the topological sort algorithm by passing an input graph and an node id map.
     * The node id mapping must be a continuous. This can be checked by setting checkMapping to
     * true.
     *
     * @param G The input graph.
     * @param nodeIdMapping Node id mapping from non-continuous to continuous ids.
     * @param checkMapping Check whether the given node id map is continuous.
     */
    GenericTopologicalSort(const GraphT &G, const NodeIdMapping &nodeIdMapping,
                           bool checkMapping = false);

    /**
     * Execute the algorithm. The algorithm is not parallel.
     */
    void run() override;

    /**
     * Return the topology
     *
     * @return One valid topology. Order in topology is from 0 to number of nodes.
     */
    const std::vector<NodeT> &getResult() const {
        assureFinished();
        return topology;
    }

private:
    enum class NodeMark : unsigned char { NONE, TEMP, PERM };

    const GraphT &G;

    std::optional<NodeIdMapping> computedNodeIdMap;

    const NodeIdMapping *nodeIdMap = nullptr;

    // Used to mark the status of each node, one vector per thread
    std::vector<NodeMark> topSortMark;

    // Contains information about the computed topology
    std::vector<NodeT> topology;

    // Helper structures
    count current;

    static NodeIdMapping computeContinuousNodeIds(const GraphT &G);

    void checkDirected();

    void checkNodeIdMap();

    index mapNode(NodeT u) const;

    // Reset algorithm data structure
    void reset();
};

using TopologicalSort = GenericTopologicalSort<Graph>;

} // namespace NetworKit

#include <networkit/graph/TopologicalSortImpl.hpp>

#endif // NETWORKIT_GRAPH_TOPOLOGICAL_SORT_HPP_
