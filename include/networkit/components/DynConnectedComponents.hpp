/*
 * DynConnectedComponents.hpp
 *
 *  Created on: June 2017
 *      Author: Eugenio Angriman
 */

#ifndef NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_

#include <map>
#include <queue>
#include <vector>

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/components/ComponentDecomposition.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup components
 * Determines and updates the connected components of an undirected graph.
 */
class DynConnectedComponents final : public ComponentDecomposition, public DynAlgorithm {

public:
    /**
     * Create ConnectedComponents class for Graph @a G.
     *
     * @param G The graph.
     */
    DynConnectedComponents(const Graph &G);

    /**
     * This method determines the connected components for the graph given in
     *  the constructor.
     */
    void run() override;

    /**
     * Updates the connected components after an edge insertion or deletion.
     *
     * @param[in]	event	The event that happened (edge deletion or
     * insertion).
     */
    void update(GraphEvent e) override;

    /**
     * Updates the connected components after a batch of edge insertions or
     * deletions.
     *
     * @param[in] batch	A vector that contains a batch of edge insertions or
     *					deletions.
     */
    void updateBatch(const std::vector<GraphEvent> &batch) override;

private:
    void addEdge(node u, node v);
    void removeEdge(node u, node v);
    void addEdgeDirected(node u, node v);
    void removeEdgeDirected(node u, node v);
    void reverseBFS(node u, node v);
    index nextAvailableComponentId(bool eraseId = true);
    void indexEdges();
    void insertEdgeIntoMap(node u, node v, edgeid eid);
    index getEdgeId(node u, node v);
    // Returns true and the corresponding edge id if the new edge was not
    // into the original graph.
    std::pair<bool, edgeid> updateMapAfterAddition(node u, node v);
    void init();
    std::pair<node, node> makePair(node u, node v);
    std::vector<bool> isTree;
    std::map<std::pair<node, node>, index> edgesMap;
    std::vector<count> tmpDistances;
    std::queue<index> componentIds;
    bool distancesInit;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
