/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMPONENTS_PARALLEL_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_PARALLEL_CONNECTED_COMPONENTS_HPP_

#include <networkit/components/ComponentDecomposition.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * @ingroup components
 * Determines the connected components of an undirected graph.
 */
class ParallelConnectedComponents final : public ComponentDecomposition {
public:
    /**
     * @param[in] G Graph for which connected components shall be computed.
     * @param[in] coarsening Specifies whether the main algorithm based on
     *   label propagation (LP) shall work recursively (true) or not (false) by
     *   coarsening/contracting an LP-computed clustering. Defaults to true
     *   since we saw positive effects in terms of running time for many networks.
     *   Beware of possible memory implications.
     */
    ParallelConnectedComponents(const Graph &G, bool coarsening = true);

    /**
     * This method determines the connected components for the graph g.
     */
    void TLX_DEPRECATED(runSequential());

    /**
     * This method determines the connected components for the graph g.
     */
    void run() override;

private:
    bool coarsening;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_PARALLEL_CONNECTED_COMPONENTS_HPP_
