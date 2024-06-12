// TODO: add boilerplate

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/SimpleHypergraphGenerator.hpp>
#include <networkit/graph/HypergraphTools.hpp>

namespace NetworKit {

SimpleHypergraphGenerator::SimpleHypergraphGenerator(count numNodes, count numEdges,
                                                     count maxEdgeOrder, bool uniform,
                                                     count regularNodeDegree)
    : numNodes{numNodes}, numEdges{numEdges}, maxEdgeOrder{maxEdgeOrder}, uniform{uniform},
      regularNodeDegree{regularNodeDegree} {

    if (regularNodeDegree != none) {
        if (uniform)
            std::runtime_error(
                "Currently, the generator does not support regular and uniform "
                "hypergraphs. Be sure only to uniform or set a regular node degree.");
        if (regularNodeDegree > numEdges)
            std::runtime_error("The node degree exceeds the number of edges. At least use as many "
                               "edges, as each node needs for it's degree.");
    }

    if (maxEdgeOrder == none)
        maxEdgeOrder = numNodes;
}

void SimpleHypergraphGenerator::satisfyEdgeOrder(Hypergraph &hGraph) {
    hGraph.forEdges([&](edgeid eid) {
        count edgeOrder = uniform == true ? maxEdgeOrder : Aux::Random::integer(maxEdgeOrder);
        hGraph.addNodesTo(HypergraphTools::randomNodes(hGraph, edgeOrder), eid);
    });
}

void SimpleHypergraphGenerator::satisfyNodeDegree(Hypergraph &hGraph) {
    hGraph.forNodes([&](node u) {
        hGraph.addNodeTo(HypergraphTools::randomEdges(hGraph, regularNodeDegree), u);
    });
}

Hypergraph SimpleHypergraphGenerator::generate() {
    // Initially the hypergraph has isolated nodes and empty hyperedges
    Hypergraph hGraph = Hypergraph(numNodes, numEdges);

    // If the hypergraph should be regular, the node degrees are used as an invariant for the
    // node-to-edge distribution
    if (regularNodeDegree != none)
        satisfyNodeDegree(hGraph);
    else
        satisfyEdgeOrder(hGraph);

    return hGraph;
}

} /* namespace NetworKit */
