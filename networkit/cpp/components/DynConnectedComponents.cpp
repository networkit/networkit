// no-networkit-format
/*
* DynConnectedComponents.cpp
*
*  Created on: June 2017
*      Author: Eugenio Angriman
*/

#include <networkit/components/DynConnectedComponents.hpp>

namespace NetworKit {

    DynConnectedComponents::DynConnectedComponents(const Graph& G) :
    ComponentDecomposition(G) {
        if (G.isDirected()) {
            throw std::runtime_error("Error, connected components of directed graphs cannot be computed, use StronglyConnectedComponents instead.");
        }
    }


    void DynConnectedComponents::init() {
        edgesMap.clear();
        component.reset(G->upperNodeIdBound(), none);
        tmpDistances.assign(G->upperNodeIdBound(), none);
        indexEdges();
        isTree.assign(edgesMap.size(), false);
        hasRun = false;
    }

    void DynConnectedComponents::run() {
        // Initializing / resetting data structures
        init();
        std::queue<node> q;
        count nComponents = 0;

        // Perform breadth-first searches
        G->forNodes([&](node u) {
            if (component[u] == none) {
                component.setUpperBound(nComponents + 1);
                index c = nComponents;
                ++nComponents;
                q.push(u);
                component[u] = c;

                do {
                    node u = q.front();
                    q.pop();
                    G->forNeighborsOf(u, [&](node v) {
                        if (component[v] == none) {
                            q.push(v);
                            component[v] = c;
                            isTree[edgesMap.at(Edge(u, v, true))] = true;
                        }
                    });
                } while (!q.empty());
            }
        });

        hasRun = true;
    }

    void DynConnectedComponents::indexEdges() {
        edgeid eid = 0;
        G->forEdges([&](node u, node v) {
            if (edgesMap.find(Edge(u, v, true)) == edgesMap.end()) {
                edgesMap.emplace(Edge{u, v, true}, eid);
                ++eid;
            }
        });
    }

    void DynConnectedComponents::update(GraphEvent event) {
        assureFinished();

        if (event.type == GraphEvent::EDGE_ADDITION) {
            addEdge(event.u, event.v);
        }
        else if (event.type == GraphEvent::EDGE_REMOVAL) {
            removeEdge(event.u, event.v);
        }
    }


    void DynConnectedComponents::updateBatch(
        const std::vector<GraphEvent>& batch
    ){
        for (auto e : batch) {
            update(e);
        }
    }


    std::pair<bool, edgeid> DynConnectedComponents::updateMapAfterAddition(
        node u, node v
    ) {
        auto it = edgesMap.find(Edge(u, v, true));

        if (it == edgesMap.end()) {
            edgeid newId = edgesMap.size();
            // Adding edge never deleted before
            edgesMap.emplace(Edge{u, v, true}, newId);
            return std::make_pair(false, none);
        }
        return std::make_pair(true, it->second);
    }


    void DynConnectedComponents::addEdge(node u, node v) {

        std::pair<bool, edgeid> updateResult = updateMapAfterAddition(u, v);

        // If u and v are already in the same component, we
        // don't have to do anything
        index minComp, maxComp;
        std::tie(minComp, maxComp) = std::minmax(component[u], component[v]);

        if (maxComp == minComp) {
            if (!updateResult.first) {
                isTree.push_back(false);
            }
            return;
        }

        // In the other case, we can merge the two components

        // All nodes in the component with higher index (maxComp) get the smaller component id (minComp).
        // The component with highest id (lastComp) is updated to maxComp, so that we can shrink
        // partition by one.
        const index lastComp = component.upperBound() - 1;
        G->parallelForNodes([&](node w) {
            if (component[w] == maxComp)
                component[w] = minComp;
            else if (component[w] == lastComp)
                component[w] = maxComp;
        });

        // lastComp is not used anymore
        component.setUpperBound(lastComp);

        if (updateResult.first) {
            isTree[updateResult.second] = true;
        }
        else {
            isTree.push_back(true);
        }
    }


    void DynConnectedComponents::removeEdge(node u, node v) {

        edgeid eid = edgesMap.at(Edge(u, v, true));

        // This edge removal does not split two components. Nothing to do.
        if (!isTree[eid]) {
            return;
        }

        isTree[eid] = false; // for coherence, we mark this edge as not valid
        std::fill(tmpDistances.begin(), tmpDistances.end(), none);

        Partition newCmp(component.getVector());
        const index nextId = newCmp.upperBound();
        newCmp.setUpperBound(nextId + 1);
        newCmp[u] = nextId;
        count newCmpSize = 0;

        std::queue<node> q;
        q.push(u);
        tmpDistances[u] = 0;

        bool connected = false;

        // Berform BFS from v to check if v reaches u
        do {
            node s = q.front();
            q.pop();

            count d = tmpDistances[s] + 1;

            // Enqueue not visited neighbors
            for (node w : G->neighborRange(s)) {
                if (tmpDistances[w] != none)
                    continue;
                tmpDistances[w] = d;
                if (w == v) { // Found another path from u to v
                    // Backtracks the path from v to u and marks all its
                    // nodes as part of the spanning tree
                    reverseBFS(u, v);
                    connected = true;
                    break;
                }

                newCmp[w] = nextId;
                ++newCmpSize;
                q.push(w);
            }
        } while (!q.empty() && !connected);

        if (!connected)
            std::swap(component, newCmp);
    }


    void DynConnectedComponents::reverseBFS(node u, node v) {

        std::queue<node> q1, q2;
        q1.push(v);

        count d = tmpDistances[v];
        count level = 1;

        do {
            do {
                const node s = q1.front();
                q1.pop();

                for (node w : G->neighborRange(s)) {
                    if (w == u) {
                        isTree[edgesMap.at(Edge(w, s, true))] = true;
                        break;
                    }

                    if ((tmpDistances[w] != none) &&
                        (d == tmpDistances[w] + level)) {
                        isTree[edgesMap.at(Edge(w, s, true))] = true;
                        q2.push(w);
                        break;
                    }
                }
            } while (!q1.empty());

            ++level;
            std::swap(q1, q2);
        } while (!q1.empty());
    }
}
