/*
 * Hypergraph.cpp
 *
 *  Created on: 24.05.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/HypergraphTools.hpp>

namespace NetworKit {

Hypergraph::Hypergraph(count n, count m, bool weighted)
    : numNodes(n), numEdges(m), nextNodeId(n), nextEdgeId(m),

      weighted(weighted), // indicates whether the graph is weighted or not
                          //   directed(directed), // indicates whether the graph is directed or not

      nodeExists(n, true), nodeWeights(weighted ? n : 0),

      nodeIncidence(n),

      edgeExists(m, true), edgeWeights(weighted ? m : 0),

      edgeIncidence(m),

      nodeAttributeMap(this), edgeAttributeMap(this) {}

node Hypergraph::addNode() {
    node v = nextNodeId; // node gets maximum id
    nextNodeId++;        // increment node range
    numNodes++;          // increment node count

    // update per node data structures
    nodeExists.push_back(true);

    nodeIncidence.emplace_back();
    if (weighted)
        nodeWeights.emplace_back();

    return v;
}

node Hypergraph::addNodes(count numberOfNewNodes) {

    nextNodeId += numberOfNewNodes;
    numNodes += numberOfNewNodes;

    // update per node data structures
    nodeExists.resize(nextNodeId, true);
    nodeIncidence.resize(nextNodeId);
    if (weighted)
        nodeWeights.resize(nextNodeId, defaultNodeWeight);

    return nextNodeId - 1;
}

node Hypergraph::addNodeTo(const std::vector<edgeid> &edges, node u) {
    if (u == none)
        u = addNode();
    for (auto eid : edges) {
        edgeIncidence[eid].insert(u);
        nodeIncidence[u].insert(eid);
    }

    return u;
}

edgeid Hypergraph::addNodesTo(const std::vector<node> &nodes, edgeid eid) {
    if (eid == none)
        eid = addEdge();
    for (auto curNode : nodes) {
        nodeIncidence[curNode].insert(eid);
        edgeIncidence[eid].insert(curNode);
    }
    return eid;
}

void Hypergraph::removeNode(node u) {
    assert(u < nextNodeId);
    assert(nodeExists[u]);

    nodeIncidence[u].clear();

    // Make the attributes of this node invalid
    auto &theMap = nodeAttributeMap.attrMap;
    for (auto it = theMap.begin(); it != theMap.end(); ++it) {
        auto attributeStorageBase = it->second.get();
        attributeStorageBase->invalidate(u);
    }

    nodeExists[u] = false;
    numNodes--;
}

void Hypergraph::restoreNode(node v) {
    assert(v < nextNodeId);
    assert(!nodeExists[v]);

    nodeExists[v] = true;
    numNodes++;
}

void Hypergraph::removeNodeFrom(node u, edgeid eid) {
    assert(eid < nextEdgeId);
    assert(edgeExists[eid]);

    edgeIncidence[eid].erase(u);
}

nodeweight Hypergraph::getNodeWeight(node u) const {
    assert(u < nextNodeId);

    nodeweight res{0.0};
    if (nodeExists[u]) {
        res = nodeWeights[u];
    }
    return res;
}

void Hypergraph::setNodeWeight(node u, nodeweight nw) {
    node tempN = u;
    if (!nodeExists[tempN])
        tempN = addNode();

    nodeWeights[tempN] = nw;
}

count Hypergraph::degree(node u) const {
    assert(u < nextNodeId);

    count res{0};
    if (nodeExists[u]) {
        res = nodeIncidence[u].size();
    }
    return res;
}

// NOTE: might profit from a parallel reduction
count Hypergraph::weightedDegree(node u) const {
    assert(u < nextNodeId);
    assert(weighted);

    count res{0};

    if (nodeExists[u]) {
        for (edgeid eid : nodeIncidence[u]) {
            res += edgeWeights[eid];
        }
    }
    return res;
}

std::set<node> Hypergraph::getNeighbors(node u) const {
    assert(u < nextNodeId);

    std::set<node> neighbors;

    for (edgeid eid : nodeIncidence[u]) {
        neighbors.insert(edgeIncidence[eid].begin(), edgeIncidence[eid].end());
    }

    neighbors.erase(u);

    return neighbors;
}

edgeid Hypergraph::addEdge() {
    edgeid eid = nextEdgeId; // edge gets maximum id
    nextEdgeId++;            // increment edge range
    numEdges++;              // increment edge count

    // update per edge data structures
    edgeExists.push_back(true);

    edgeIncidence.emplace_back();
    if (weighted)
        edgeWeights.emplace_back();

    return eid;
}

edgeid Hypergraph::addEdge(const std::vector<node> &nodes, bool addMissing) {

    edgeid eid = addEdge();
    edgeIncidence[eid] = std::set<node>(nodes.begin(), nodes.end());

    if (addMissing) {
        for (auto v : edgeIncidence[eid]) {
            node currentMax = numNodes;
            while (v > currentMax) {
                currentMax = addNode();
                nodeExists[currentMax] = false;
            }
            addNode();
        }
    }

    for (auto v : edgeIncidence[eid]) {
        nodeIncidence[v].insert(eid);
    }

    return eid;
}

void Hypergraph::removeEdge(edgeid eid) {
    assert(eid < nextEdgeId);
    assert(edgeExists[eid]);

    edgeIncidence[eid].clear();

    // Make the attributes of this edge invalid
    auto &theMap = edgeAttributeMap.attrMap;
    for (auto it = theMap.begin(); it != theMap.end(); ++it) {
        auto attributeStorageBase = it->second.get();
        attributeStorageBase->invalidate(eid);
    }

    edgeExists[eid] = false;
    numEdges--;
}

edgeweight Hypergraph::getEdgeWeight(edgeid eid) const {
    assert(eid < nextEdgeId);
    return weighted ? edgeWeights[eid] : defaultEdgeWeight;
}

void Hypergraph::setEdgeWeight(edgeid eid, edgeweight ew) {
    if (!weighted) {
        throw std::runtime_error("Cannot set edge weight in unweighted hypergraphs.");
    }
    edgeid tempEid = eid;
    if (!edgeExists[tempEid])
        tempEid = addEdge();

    edgeWeights[tempEid] = ew;
}
} // namespace NetworKit
