// TODO: add boilerplate

#include <ranges>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/HypergraphTools.hpp>

namespace NetworKit {

node HypergraphTools::randomNode(const Hypergraph &hGraph) {
    if (!hGraph.numberOfNodes())
        return none;

    auto &gen = Aux::Random::getURNG();
    std::uniform_int_distribution<node> distr{0, hGraph.upperNodeIdBound() - 1};
    node v;

    do {
        // When there are many deleted nodes, we might call Aux::Random::integer
        // many times, and it is very expensive.
        v = distr(gen);
    } while (!hGraph.hasNode(v));

    return v;
}

std::vector<node> HypergraphTools::randomNodes(const Hypergraph &hGraph, count numNodes) {
    assert(numNodes <= hGraph.numberOfNodes());
    std::vector<node> selectedNodes;
    std::vector<bool> alreadySelected(hGraph.numberOfNodes(), false);

    if (numNodes == hGraph.numberOfNodes()) {
        selectedNodes.insert(selectedNodes.begin(), hGraph.nodeRange().begin(),
                             hGraph.nodeRange().end());
    } else if (numNodes
               > hGraph.numberOfNodes() / 2) { // in order to minimize the calls to randomNode
                                               // we randomize the ones that aren't pivot
                                               // if the are more to be selected than not-selected
        for (count i = 0; i < hGraph.numberOfNodes() - numNodes;
             ++i) { // we have to sample distinct nodes
            node v = HypergraphTools::randomNode(hGraph);
            while (alreadySelected[v]) {
                v = HypergraphTools::randomNode(hGraph);
            }
            alreadySelected[v] = true;
        }

        for (const auto sample : hGraph.nodeRange()) {
            if (!alreadySelected[sample]) { // recall that we selected the non-pivot nodes
                selectedNodes.push_back(sample);
                if (selectedNodes.size() == numNodes)
                    break;
            }
        }
    } else {
        for (count i = 0; i < numNodes; ++i) { // we have to selected distinct nodes
            node v = HypergraphTools::randomNode(hGraph);
            while (alreadySelected[v]) {
                v = HypergraphTools::randomNode(hGraph);
            }
            selectedNodes.push_back(v);
            alreadySelected[v] = true;
        }
    }
    return selectedNodes;
}

edgeid HypergraphTools::randomEdge(const Hypergraph &hGraph) {
    if (!hGraph.numberOfEdges())
        return none;

    auto &gen = Aux::Random::getURNG();
    std::uniform_int_distribution<edgeid> distr{0, hGraph.upperEdgeIdBound() - 1};
    edgeid eid;

    do {
        // When there are many deleted nodes, we might call Aux::Random::integer
        // many times, and it is very expensive.
        eid = distr(gen);
    } while (!hGraph.hasEdge(eid));

    return eid;
}

std::vector<edgeid> HypergraphTools::randomEdges(const Hypergraph &hGraph, count numEdges) {
    assert(numEdges <= hGraph.numberOfEdges());
    std::vector<edgeid> selectedEdges;
    std::vector<bool> alreadySelected(hGraph.numberOfEdges(), false);

    if (numEdges == hGraph.numberOfEdges()) {
        selectedEdges.resize(hGraph.numberOfEdges(), 0);
        index i = 0;
        hGraph.forEdges([&](edgeid eId) {
            if (hGraph.hasEdge(eId)) {
                selectedEdges[i] = eId;
                i++;
            }
        });
    } else if (numEdges
               > hGraph.numberOfEdges() / 2) { // in order to minimize the calls to randomEdge
                                               // we randomize the ones that aren't pivot
                                               // if the are more to be selected than not-selected
        for (count i = 0; i < hGraph.numberOfEdges() - numEdges;
             ++i) { // we have to sample distinct edges
            edgeid sample = HypergraphTools::randomEdge(hGraph);
            while (alreadySelected[sample]) {
                sample = HypergraphTools::randomEdge(hGraph);
            }
            selectedEdges.push_back(sample);
            alreadySelected[sample] = true;
        }

        for (edgeid eId = 0; eId < hGraph.numberOfEdges(); ++eId) {
            if (hGraph.hasEdge(eId) && !alreadySelected[eId]) {
                selectedEdges.push_back(eId);
                INFO("Reached here: ", selectedEdges.size());
                if (selectedEdges.size() == numEdges) {
                    INFO("Reached here2: ", selectedEdges.size());
                    break;
                }
            }
        };
    } else {
        for (count i = 0; i < numEdges; ++i) { // we have to selected distinct nodes
            edgeid sample = HypergraphTools::randomEdge(hGraph);
            while (alreadySelected[sample]) {
                sample = HypergraphTools::randomEdge(hGraph);
            }
            selectedEdges.push_back(sample);
            alreadySelected[sample] = true;
        }
    }
    return selectedEdges;
}

count HypergraphTools::maxEdgeOrder(const Hypergraph &hGraph) {
    count result = 0;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(max : result)
    for (omp_index eid = 0; eid < static_cast<omp_index>(hGraph.upperEdgeIdBound()); ++eid) {
        result = std::max(result, hGraph.order(static_cast<edgeid>(eid)));
    }
#else
    hGraph.forEdges([&](edgeid eid) { result = std::max(result, hGraph.order(eid)); });
#endif
    return result;
}

count HypergraphTools::maxDegree(const Hypergraph &hGraph) {
    count result = 0;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(max : result)
    for (omp_index u = 0; u < static_cast<omp_index>(hGraph.upperNodeIdBound()); ++u) {
        result = std::max(result, hGraph.degree(static_cast<node>(u)));
    }
#else
    hGraph.forNodes([&](node u) { result = std::max(result, hGraph.degree(u)); });
#endif
    return result;
}

edgeweight HypergraphTools::maxWeightedDegree(const Hypergraph &hGraph) {
    edgeweight result = 0;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(max : result)
    for (omp_index u = 0; u < static_cast<omp_index>(hGraph.upperNodeIdBound()); ++u) {
        result = std::max(result, hGraph.weightedDegree(static_cast<node>(u)));
    }
#else
    hGraph.forNodes([&](node u) { result = std::max(result, hGraph.weightedDegree(u)); });
#endif
    return result;
}

std::unordered_set<node> HypergraphTools::getIntersection(Hypergraph &hGraph, edgeid eid1,
                                                          edgeid eid2) {
    std::unordered_set<node> smallerUSet = hGraph.order(eid1) < hGraph.order(eid2)
                                               ? hGraph.edgeMembers(eid1)
                                               : hGraph.edgeMembers(eid2);
    std::unordered_set<node> largerUSet = hGraph.order(eid1) < hGraph.order(eid2)
                                              ? hGraph.edgeMembers(eid2)
                                              : hGraph.edgeMembers(eid1);

    std::unordered_set<node> uSetIntersection;

    for (node u : smallerUSet) {
        if (largerUSet.count(u) > 0)
            uSetIntersection.insert(u);
    }

    return uSetIntersection;
}

count HypergraphTools::getIntersectionSize(Hypergraph &hGraph, edgeid eid1, edgeid eid2) {
    std::unordered_set<node> smallerUSet = hGraph.order(eid1) < hGraph.order(eid2)
                                               ? hGraph.edgeMembers(eid1)
                                               : hGraph.edgeMembers(eid2);
    std::unordered_set<node> largerUSet = hGraph.order(eid1) < hGraph.order(eid2)
                                              ? hGraph.edgeMembers(eid2)
                                              : hGraph.edgeMembers(eid1);

    count intersectionSize = 0;

    for (node u : smallerUSet) {
        if (largerUSet.count(u) > 0)
            intersectionSize++;
    }

    return intersectionSize;
}

Graph HypergraphTools::cliqueExpansion(Hypergraph &hGraph) {

    Graph cliqueExpansion(hGraph.numberOfNodes());

    hGraph.forEdges([&](edgeid eid) {
        const std::unordered_set<node> &nodesInEdge = hGraph.edgeMembers(eid);
        for (auto firstIt = nodesInEdge.begin(); firstIt != nodesInEdge.end(); ++firstIt) {
            for (auto secondIt = std::next(firstIt); secondIt != nodesInEdge.end(); ++secondIt) {
                node v = *firstIt;
                node w = *secondIt;
                if (v > w)
                    cliqueExpansion.addEdge(v, w);
            }
        }
    });

    cliqueExpansion.removeMultiEdges();

    return cliqueExpansion;
}

Graph HypergraphTools::lineExpansion(Hypergraph &hGraph) {
    std::map<std::pair<node, edgeid>, node> lineMap;

    count expansionSize = 0;
    hGraph.forEdges([&](edgeid eid) { expansionSize += hGraph.order(eid); });

    // now create the lineExpansion graph since we now know the amount of nodes in it
    Graph lineExpansion(expansionSize);
    auto attrNodeRef = lineExpansion.nodeAttributes().attach<node>("node");
    auto attrEdgeRef = lineExpansion.nodeAttributes().attach<edgeid>("edgeid");

    // First add all edges inside each clique
    // In addition set the node attributes in order to maintain the original data
    count currentId = 0;
    hGraph.forEdges([&](edgeid eid) {
        const std::unordered_set<node> &nodesInEdge = hGraph.edgeMembers(eid);
        count currentOrder = hGraph.order(eid);
        node offset = 0;
        for (auto firstIt = nodesInEdge.begin(); firstIt != nodesInEdge.end(); ++firstIt) {
            attrNodeRef[currentId + offset] = node{*firstIt};
            attrEdgeRef[currentId + offset] = edgeid{eid};
            lineMap[std::make_pair(node{*firstIt} + offset, edgeid{eid})] = currentId + offset;
            for (count nextIdsInEdge = currentId + offset; nextIdsInEdge < currentId + currentOrder;
                 ++nextIdsInEdge) {
                lineExpansion.addEdge(currentId, nextIdsInEdge);
            }
            offset++;
        }
        currentId += currentOrder;
    });

    // Now add all edges between hyperedges
    hGraph.forNodes([&](node u) {
        const std::unordered_set<node> &edgesOfNode = hGraph.edgesOf(u);
        for (auto firstIt = edgesOfNode.begin(); firstIt != edgesOfNode.end(); ++firstIt) {
            for (auto secondIt = std::next(firstIt); secondIt != edgesOfNode.end(); ++secondIt) {
                node v = lineMap[std::make_pair(u, *firstIt)];
                node w = lineMap[std::make_pair(u, *secondIt)];
                lineExpansion.addEdge(v, w);
            }
        }
    });

    lineExpansion.removeMultiEdges();

    return lineExpansion;
}

    return std::unordered_set<node>(view.begin(), view.end());
}

} // namespace NetworKit
