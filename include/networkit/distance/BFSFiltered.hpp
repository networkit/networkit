/*
 * BFSFiltered.hpp
 *
 *  Created on: Jan 21, 2024
 *      Author: Felix
 */

#ifndef NETWORKIT_DISTANCE_BFSFILTERED_HPP_
#define NETWORKIT_DISTANCE_BFSFILTERED_HPP_

#include <vector>

#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * The BFSFiltered class is used to do a breadth-first search on a Graph from a given
 * source node.
 */
class BFSFiltered final : public SSSP {

public:
    /**
     * Constructs the BFSFiltered class for @a G and source node @a source.
     *
     * @param G The graph
     * @param source The source node of the breadth-first search
     * @param storePaths Paths are reconstructable and the number of paths is
     * stored.
     * @param filteredFirst, filteredLast Range of nodes to be ignored
     * @param storeNodesSortedByDistance Store a vector of nodes ordered in
     * increasing distance from the source.
     * @param target The target node.
     */
    template <class InputIt>
    BFSFiltered(const Graph &G, node source, InputIt filteredFirst, InputIt filteredLast, bool storePaths = true,
        bool storeNodesSortedByDistance = false, node target = none): SSSP(G, source, storePaths, storeNodesSortedByDistance, target), filteredNodes(filteredFirst, filteredLast) {}

    /**
     * Breadth-first search from @a source.
     */
    void run() override;

protected:
    const std::vector<node> filteredNodes;
};
} /* namespace NetworKit */
#endif // NETWORKIT_DISTANCE_BFSFILTERED_HPP_
