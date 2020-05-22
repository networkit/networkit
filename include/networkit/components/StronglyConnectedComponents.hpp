/*
 * StronglyConnectedComponents.hpp
 *
 *  Created on: 01.06.2014
 *      Authors: Klara Reichard <klara.reichard@gmail.com>
 *               Marvin Ritter <marvin.ritter@gmail.com>
 *               Obada Mahdi <omahdi@gmail.com>
 *               Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_COMPONENTS_STRONGLY_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_STRONGLY_CONNECTED_COMPONENTS_HPP_

#include <map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * @ingroup components
 *
 */
class StronglyConnectedComponents final : public Algorithm {

public:
    /**
     * Determines the strongly connected components of a directed graph using Tarjan's algorithm.
     *
     * @param G Graph A directed graph.
     */
    StronglyConnectedComponents(const Graph &G);

    /**
     * Determines the strongly connected components of a directed graph using Tarjan's algorithm.
     *
     * @param G Graph A directed graph.
     * @param iterativeAlgo bool Whether to use the iterative algorithm or the recursive one.
     *
     * This constructor has been deprecated because the recursive algorithm has been deprecated.
     */
    TLX_DEPRECATED(StronglyConnectedComponents(const Graph &G, bool iterativeAlgo));

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * This method determines the connected components for the graph g
     * (iterative implementation).
     *
     * This method is deprecated, run() already uses the iterative implementation.
     */
    void TLX_DEPRECATED(runIteratively());

    /**
     * This method determines the connected components for the graph g
     * (recursive implementation).
     *
     * This method is deprecated, use run().
     */
    void TLX_DEPRECATED(runRecursively());

    /**
     * Return the number of connected components.
     *
     * @return count The number of components of the graph.
     */
    count numberOfComponents() const {
        assureFinished();
        return componentIndex;
    }

    /**
     * Returns the component of the input node.
     *
     * @param[in] u The input node.
     * @return count The index of the component of the input node.
     */
    index componentOfNode(node u) const {
        assureFinished();
        assert(G->hasNode(u));
        return component[u];
    }

    /**
     * Return a Partition object that represents the components.
     *
     * @return Partition Partitioning of the strongly connected components.
     */
    Partition getPartition() const {
        assureFinished();
        return Partition(component);
    }

    /**
     * Return a map with the component indexes as keys, and their size as values.
     *
     * @return std::map<index, count> Map with components indexes as keys, and their size as values.
     */
    std::map<node, count> getComponentSizes() const {
        assureFinished();

        std::vector<count> componentSizes(componentIndex, 0);
        G->forNodes([&](const node u) { ++componentSizes[component[u]]; });

        std::map<node, count> result;
        for (index i = 0; i < componentIndex; ++i)
            result[i] = componentSizes[i];

        return result;
    }

    /**
     * Return a list of components.
     *
     * @return std::vector<std::vector<node>> List of components.
     */
    std::vector<std::vector<node>> getComponents() const {
        assureFinished();
        std::vector<std::vector<node>> result(componentIndex);

        G->forNodes([&](const node u) { result[component[u]].push_back(u); });

        return result;
    }

private:
    const Graph *G;
    bool iterativeAlgo;
    std::vector<index> component;
    index componentIndex;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_STRONGLY_CONNECTED_COMPONENTS_HPP_
