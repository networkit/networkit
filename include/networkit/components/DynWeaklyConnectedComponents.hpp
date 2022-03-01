/*
 * DynWeaklyConnectedComponents.hpp
 *
 *  Created on: June 20, 2017
 *      Author: Eugenio Angriman
 */

#ifndef NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_

#include <unordered_map>
#include <vector>

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/components/ComponentDecomposition.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup components
 * Determines and updates the weakly connected components of a directed
 * graph.
 */
class DynWeaklyConnectedComponents final : public ComponentDecomposition, public DynAlgorithm {

public:
    /**
     * Create DynWeaklyConnectedComponents class for Graph @a G.
     *
     * @param G The graph.
     */
    DynWeaklyConnectedComponents(const Graph &G);

    /**
     * This method determines the weakly connected components for the graph
     * given in the constructor.
     */
    void run() override;

    /**
     * Updates the weakly connected components after an edge insertion or
     * deletion.
     *
     * @param[in]	event	The event that happened (edge insertion
     *                       or deletion).
     */
    void update(GraphEvent event) override;

    /**
     * Updates the weakly connected components after a batch of edge events
     * (insertions or deletions).
     *
     * @param[in] batch  A vector that contains a batch of edge events
     *                   (insertions or deletions).
     */
    void updateBatch(const std::vector<GraphEvent> &batch) override;

private:
    void init();
    void addEdge(node u, node v);
    void removeEdge(node u, node v);
    void indexEdges();
    edgeid updateMapAfterAddition(node u, node v);
    void updateTreeAfterAddition(edgeid eid, bool partOfTree);
    void reverseBFS(node u, node v);

    // Whether each edge is part of a spanning tree. If an edge is removed
    // then a components could split.
    std::vector<bool> isTree;

    // Maps the edges to their IDs
    std::unordered_map<Edge, edgeid> edgesMap;

    // Temporary distances used for BFSs after edge updates.
    std::vector<count> tmpDistances;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_
