/*
 * DynConnectedComponents.hpp
 *
 *  Created on: June 2017
 *      Author: Eugenio Angriman
 */

#ifndef NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_

#include <memory>
#include <vector>

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

    ~DynConnectedComponents() override;

    /**
     * Finds the connected components of the input graph.
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
    std::unique_ptr<DynConnectedComponentsDetails::DynConnectedComponentsImpl<false>> impl;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
