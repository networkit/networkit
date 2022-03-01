// no-networkit-format
/*
* DynDynWeaklyConnectedComponents.cpp
*
*  Created on: June 20, 2017
*      Author: Eugenio Angriman
*/

#include <networkit/components/DynWeaklyConnectedComponents.hpp>

namespace NetworKit {

    DynWeaklyConnectedComponents::DynWeaklyConnectedComponents(const Graph& G) : ComponentDecomposition(G) {
        if (!G.isDirected()) {
            throw std::runtime_error("Weakly Connected Components can be computeed for directed graphs. Use ConnectedComponents for undirected graphs.");
        }
    }


    void DynWeaklyConnectedComponents::init() {
        edgesMap.clear();
        indexEdges();
        component.reset(G->upperNodeIdBound(), none);
        isTree.assign(edgesMap.size(), false);
        hasRun = false;
    }


    void DynWeaklyConnectedComponents::run() {

        init();

        // Queue for BFS.
        std::queue<node> q;

        count nComponents = 0;

        // Perform BFSs to assign a component ID to each node.
        G->forNodes([&](node u) {

            // Node u has not been visited.
            if (component[u] == none) {

                component.setUpperBound(nComponents + 1);
                index c = nComponents;
                ++nComponents;
                component[u] = c;

                // Start a new BFS from u.
                q.push(u);

                do {
                    node v = q.front();
                    q.pop();

        
                    auto updateComponent = [&](node, node u, edgeweight, edgeid) -> void {
                        if (component[u] == none) {
                            q.push(u);
                            component[u] = c;
                            isTree[edgesMap.at(Edge(u, v, true))] = true;
                        }
                    };

                    // Enqueue neighbors (both from in and out edges) and set
                    // new component.
                    G->forNeighborsOf(v, updateComponent);
                    G->forInNeighborsOf(v, updateComponent);
                } while (!q.empty());
            }
        });

        hasRun = true;
    }

    void DynWeaklyConnectedComponents::indexEdges() {
        edgeid eid = 0;
        G->forEdges([&] (node u, node v) {
            if (edgesMap.find(Edge(u, v, true)) == edgesMap.end()) {
                edgesMap.emplace(Edge{u, v, true}, eid);
                ++eid;
            }
        });
    }

    edgeid DynWeaklyConnectedComponents::updateMapAfterAddition(
        node u, node v
    ) {

        Edge newEdge(u, v, true);
        auto it = edgesMap.find(newEdge);

        if (it == edgesMap.end()) {
            // Adding edge never deleted before
            index newId = edgesMap.size();
            edgesMap.emplace(newEdge, newId);
            return newId;
        }

        return it->second;
    }


    void DynWeaklyConnectedComponents::update(GraphEvent event) {

        assureFinished();

        if (event.type == GraphEvent::EDGE_ADDITION) {
            addEdge(event.u, event.v);
        }
        else if (event.type == GraphEvent::EDGE_REMOVAL) {
            removeEdge(event.u, event.v);
        }
    }


    void DynWeaklyConnectedComponents::updateBatch(
        const std::vector<GraphEvent> &
    ) {
        run();
    }


    void DynWeaklyConnectedComponents::addEdge(node u, node v) {

        edgeid eid = updateMapAfterAddition(u, v);

        // If u and v are already in the same component, we
        // don't have to do anything. Same thing if edge (v, u) already exists.
        index minComp, maxComp;
        std::tie(minComp, maxComp) = std::minmax(component[u], component[v]);

        if (maxComp == minComp || G->hasEdge(v, u)) {
            updateTreeAfterAddition(eid, false);
            return;
        }

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
        updateTreeAfterAddition(eid, true);
    }


    void DynWeaklyConnectedComponents::updateTreeAfterAddition(
        edgeid eid, bool partOfTree
    ) {
        if (isTree.size() > eid) {
            isTree[eid] = partOfTree;
        }
        else if (isTree.size() == eid) {
            isTree.push_back(partOfTree);
        }
        else throw std::runtime_error("Edge indexing error");
    }


    void DynWeaklyConnectedComponents::removeEdge(node u, node v) {

        edgeid eid = edgesMap.at(Edge(u, v, true));

        // If (u, v) is not part of the spanning tree or if edge (v, u) already
        // exists we don't have to do nothing.
        if (!isTree[eid]) {
            return;
        }

        // Edge "eid" is removed from the graph. For performance reasons we
        // keep it in memory, for coherence we claim that it is no more part of
        // the spanning tree.
        isTree[eid] = false;

        Partition newCmp(component.getVector());
        const index nextId = newCmp.upperBound();
        newCmp.setUpperBound(nextId + 1);
        newCmp[u] = nextId;
        count newCmpSize = 1;

        tmpDistances.assign(G->upperNodeIdBound(), none);
        std::queue<node> q;
        q.push(u);
        tmpDistances[u] = 0;

        bool connected = false;

        // Berform BFS from v to check if v reaches u
        do {
            node s = q.front();
            q.pop();

            count d = tmpDistances[s] + 1;

            // Enqueue not visited neighbors reachable through outgoing edges
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

            // Checking in neighbors (forNeighborsOf gets only the out-degree
            // neighbors).
            if (connected)
                break;

            for (node w : G->inNeighborRange(s)) {
                if (tmpDistances[w] != none)
                    continue;
                tmpDistances[w] = d;
                if (w == v) { // Found another path from u to v
                    // Backtracks the path from v to u and marks all
                    // its nodes as part of the spanning tree
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

    void DynWeaklyConnectedComponents::reverseBFS(node u, node v) {

        std::queue<node> q1, q2;
        q1.push(v);
        count d = tmpDistances[v];
        count level = 1;

        do {
            do {
                const node s = q1.front();
                q1.pop();

                bool nextEdgeFound = false;

                auto visitNeighborsReversed = [&](node w) -> bool {
                    // Reverse BFS finished
                    if (w == u) {
                        isTree[edgesMap.at(Edge(w, s, true))] = true;
                        nextEdgeFound = true;
                        return true;
                    }

                    // Found next node for reverse BFS
                    if ((tmpDistances[w] != none) && (d == tmpDistances[w] + level)) {
                        isTree[edgesMap.at(Edge(w, s, true))] = true;
                        nextEdgeFound = true;
                        q2.push(w);
                        return true;
                    }

                    // Discarding node from reverse path
                    return false;
                };

                for (node w : G->neighborRange(s))
                    if (visitNeighborsReversed(w))
                        break;

                if (nextEdgeFound)
                    break;

                for (node w : G->inNeighborRange(s))
                    if (visitNeighborsReversed(w))
                        break;
            } while (!q1.empty());
            std::swap(q1, q2);
            ++level;
        } while (!q1.empty());
    }
}
