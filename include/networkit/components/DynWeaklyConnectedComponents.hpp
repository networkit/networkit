/*
 * DynWeaklyConnectedComponents.hpp
 *
 *  Created on: June 20, 2017
 *      Author: Eugenio Angriman
 */

#ifndef NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_

#include <memory>

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/components/ComponentDecomposition.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

// pImpl
namespace DynConnectedComponentsDetails {
template <bool>
class DynConnectedComponentsImpl;
} // namespace DynConnectedComponentsDetails

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

    ~DynWeaklyConnectedComponents() override;

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
    std::unique_ptr<DynConnectedComponentsDetails::DynConnectedComponentsImpl<true>> impl;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_
