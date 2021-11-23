/*
 * GraphCoarsening.hpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_COARSENING_GRAPH_COARSENING_HPP_
#define NETWORKIT_COARSENING_GRAPH_COARSENING_HPP_

#include <map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup coarsening
 * Abstract base class for graph coarsening/contraction algorithms.
 */
class GraphCoarsening : public Algorithm {

public:
    GraphCoarsening(const Graph &G);

    ~GraphCoarsening() override = default;

    void run() override = 0;

    const Graph &getCoarseGraph() const;

    Graph &getCoarseGraph();

    /**
     * Get mapping from fine to coarse node.
     */
    const std::vector<node> &getFineToCoarseNodeMapping() const;

    std::vector<node> &getFineToCoarseNodeMapping();

    /**
     * Get mapping from coarse node to collection of fine nodes.
     */
    std::map<node, std::vector<node>> getCoarseToFineNodeMapping() const;

protected:
    const Graph *G;
    Graph Gcoarsened;
    std::vector<node> nodeMapping;
};

} // namespace NetworKit

#endif // NETWORKIT_COARSENING_GRAPH_COARSENING_HPP_
