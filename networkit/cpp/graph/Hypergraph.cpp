/*
 * Hypergraph.cpp
 *
 *  Created on: 24.05.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

Hypergraph::Hypergraph(count n, count m, bool weighted, bool directed)
    : n(n), m(m), z(n), omega(0),

      weighted(weighted), // indicates whether the graph is weighted or not
                          //   directed(directed), // indicates whether the graph is directed or not

      nodeExists(n, true), nodeWeights(weighted ? n : 0),

      nodeIncidence(directed ? n : 0),

      edgeExists(m, true), edgeWeights(weighted ? m : 0),

      edgeIncidence(directed ? m : 0),

      nodeAttributeMap(this), edgeAttributeMap(this) {}

node Hypergraph::addNode() {
    node v = z; // node gets maximum id
    z++;        // increment node range
    n++;        // increment node count

    // update per node data structures
    nodeExists.push_back(true);

    nodeIncidence.emplace_back();
    if (weighted)
        nodeWeights.emplace_back();
}

return v;
}

node Hypergraph::addNodeTo(std::vector<edgeid> edges, node u) {
    if (u == none)
        u = addNode();
    for (auto eid : edges) {
        edgeIncidence[eid].insert(u);
    }

    return v;
}

void Hypergraph::removeNode(node v) {
    assert(v < z);
    assert(nodeExists[v]);

    nodeIncidence.clear();

    // Make the attributes of this node invalid
    auto &theMap = nodeAttributeMap.attrMap;
    for (auto it = theMap.begin(); it != theMap.end(); ++it) {
        auto attributeStorageBase = it->second.get();
        attributeStorageBase->invalidate(v);
    }

    exists[v] = false;
    n--;
}

void Hypergraph::removeNodeFrom(node v, edgeid eid) {
    assert(eid < omega);
    assert(edgeExists[eid]);

    edgeIncidence[eid].remove(v);
}

edgeid Hypergraph::addEdge() {
    edgeid eid = omega; // edge gets maximum id
    omega++;            // increment edge range
    m++;                // increment edge count

    // update per edge data structures
    edgeExists.push_back(true);

    edgeIncidence.emplace_back();
    if (weighted)
        edgeWeights.emplace_back();
}
}

} // namespace NetworKit
