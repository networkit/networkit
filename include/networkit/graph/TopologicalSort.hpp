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
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Given a directed graph G, the topology sort algorithm creates one valid topology order of nodes.
 * Undirected graphs are not accepted as input, since a topology sort is a linear ordering of
 * vertices such that for every edge u -> v, node u comes before v in the ordering.
 */
class TopologicalSort final : public Algorithm {
public:
    /**
     * Initialize the topological sort algorithm by passing an input graph. Note that topological
     * sort is defined for directed graphs only. Node ids must be continuous in the interval [0, n).
     *
     * @param G The input graph.
     */
    TopologicalSort(const Graph &G);

    /**
     * Initialize the topological sort algorithm by passing an input graph. Note that topological
     * sort is defined for directed graphs only. The node id mapping must be a continuous.
     *
     * @param G The input graph.
     */
    TopologicalSort(const Graph &G, std::unordered_map<node, node> &nodeIdMapping,
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
    const std::vector<node> &getResult() const {
        assureFinished();
        return topology;
    }

private:
    enum class NodeMark : unsigned char { NONE, TEMP, PERM };

    const Graph *G;

    std::optional<std::unordered_map<node, node>> computedNodeIdMap;
    
    std::unordered_map<node, node> *nodeIdMap;

    // Used to mark the status of each node, one vector per thread
    std::vector<NodeMark> topSortMark;

    // Contains information about the computed topology
    std::vector<node> topology;

    // Helper structures
    count current;

    void checkDirected();

    // Reset algorithm data structure
    void reset();
};
} // namespace NetworKit

#endif // NETWORKIT_GRAPH_TOPOLOGICAL_SORT_HPP_
