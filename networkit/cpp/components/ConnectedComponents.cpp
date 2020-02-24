/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

// networkit-format

#include <set>
#include <unordered_map>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

ConnectedComponents::ConnectedComponents(const Graph &G) : G(&G) {
    if (G.isDirected()) {
        throw std::runtime_error("Error, connected components of directed graphs cannot be "
                                 "computed, use StronglyConnectedComponents for them.");
    }
}

void ConnectedComponents::run() {
    DEBUG("initializing labels");
    component = Partition(G->upperNodeIdBound(), none);
    numComponents = 0;

    std::queue<node> q;

    // perform breadth-first searches
    G->forNodes([&](node u) {
        if (component[u] == none) {
            component.setUpperBound(numComponents + 1);
            index c = numComponents;

            q.push(u);
            component[u] = c;

            do {
                node u = q.front();
                q.pop();
                // enqueue neighbors, set component
                G->forNeighborsOf(u, [&](node v) {
                    if (component[v] == none) {
                        q.push(v);
                        component[v] = c;
                    }
                });
            } while (!q.empty());

            ++numComponents;
        }
    });

    hasRun = true;
}

Partition ConnectedComponents::getPartition() const {
    assureFinished();
    return this->component;
}

std::vector<std::vector<node>> ConnectedComponents::getComponents() const {
    assureFinished();

    // transform partition into vector of unordered_set
    std::vector<std::vector<node>> result(numComponents);

    G->forNodes([&](node u) { result[component[u]].push_back(u); });

    return result;
}

std::map<index, count> ConnectedComponents::getComponentSizes() const {
    assureFinished();
    return this->component.subsetSizeMap();
}

Graph ConnectedComponents::extractLargestConnectedComponent(const Graph &G, bool compactGraph) {
    if (!G.numberOfNodes()) {
        return G;
    }

    ConnectedComponents cc(G);
    cc.run();

    const auto compSizes = cc.getComponentSizes();
    if (compSizes.size() == 1) {
        if (compactGraph) {
            return GraphTools::getCompactedGraph(G, GraphTools::getContinuousNodeIds(G));
        }
        return G;
    }

    const auto largestCC =
        std::max_element(compSizes.begin(), compSizes.end(),
                         [](const std::pair<index, count> &x, const std::pair<index, count> &y) {
                             return x.second < y.second;
                         });

    if (compactGraph) {
        std::unordered_map<node, node> continuousNodeIds;
        index nextId = 0;
        G.forNodes([&](const node u) {
            if (cc.componentOfNode(u) == largestCC->first) {
                continuousNodeIds[u] = nextId++;
            }
        });

        return GraphTools::getRemappedGraph(
            G, largestCC->second, [&](const node u) { return continuousNodeIds[u]; },
            [&](const node u) { return cc.componentOfNode(u) != largestCC->first; });

    } else {
        Graph S(G);
        const auto components = cc.getComponents();
        for (size_t i = 0; i < components.size(); ++i) {
            if (i != largestCC->first) {
                for (const auto u : components[i]) {
                    S.removeNode(u);
                }
            }
        }

        return S;
    }
}

} // namespace NetworKit
