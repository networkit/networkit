/*
 * SimpleHypergraphGenerator.hpp
 *
 *  Created on: 12.06.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#ifndef NETWORKIT_GENERATORS_SIMPLEHYPERGRAPHGENERATOR_HPP_
#define NETWORKIT_GENERATORS_SIMPLEHYPERGRAPHGENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * A simple hypergraph generator, unweighted, uniform and regular hypergraphs.
 */
class SimpleHypergraphGenerator final : public StaticGraphGenerator<Hypergraph> {

public:
    /**
     * Creates simple random hyper graphs.
     *
     * The generator can be used to generate unweighted hypergraphs with a randomly
     * distributed edge order and node degree. The maximum order/degree can be controlled by
     * setting parameters maxEdgeOrder/maxNodeDegree. Optionally, these bounds can be
     * set to be the same for every node (regular) or edge (uniform).
     *
     * @param numNodes  Number of nodes n in the hypergraph.
     * @param numEdges  Number of nodes n in the hypergraph.
     * @param maxEdgeOrder  Uppper bound for the edge order.
     * @param uniform   Controls, whether the hypergraph is uniform.
     * @param regularNodeDegree Constant node degree.
     */
    SimpleHypergraphGenerator(count numNodes, count numEdges, count maxEdgeOrder = none,
                              bool uniform = false, count regularNodeDegree = none);

    Hypergraph generate() override;

private:
    void satisfyEdgeOrder(Hypergraph &hGraph);

    void satisfyNodeDegree(Hypergraph &hGraph);

    count numNodes;
    count numEdges;
    count maxEdgeOrder;
    bool uniform;
    count regularNodeDegree;
};

} // namespace NetworKit

#endif //
