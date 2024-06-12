// TODO: add boilerplate

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

} // namespace NetworKit
