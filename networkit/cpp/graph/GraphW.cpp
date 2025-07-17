/*
 * GraphW.cpp
 *
 *  Created on: 2025-07-16
 *      Author: Christian Staudt
 *              Klara Reichard <klara.reichard@gmail.com>
 *              Marvin Ritter <marvin.ritter@gmail.com>
 *              Arun Sharma
 */

#include <networkit/graph/GraphW.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <sstream>
#include <unordered_set>
#include <ranges>

namespace NetworKit {

GraphW::GraphW(std::initializer_list<WeightedEdge> edges) : Graph(0, true) {
    using namespace std;

    /* Number of nodes = highest node index + 1 */
    for (const auto &edge : edges) {
        node x = std::max(edge.u, edge.v);
        while (x >= z) {
            addNode();
        }
    }

    /* Add edges */
    for (const auto &edge : edges) {
        addEdge(edge.u, edge.v, edge.weight);
    }
}

/** EDGE IDS **/

void GraphW::indexEdges(bool force) {
    if (edgesIndexed && !force)
        return;

    omega = 0; // reset edge ids (for re-indexing)

    outEdgeIds.clear(); // reset ids vector (for re-indexing)
    outEdgeIds.resize(outEdges.size());
    forNodes([&](node u) { outEdgeIds[u].resize(outEdges[u].size(), none); });

    if (directed) {
        inEdgeIds.resize(inEdges.size());
        forNodes([&](node u) { inEdgeIds[u].resize(inEdges[u].size(), none); });
    }

    // assign edge ids for edges in one direction
    forNodes([&](node u) {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            node v = outEdges[u][i];
            if (v != none && (directed || (u >= v))) {
                // new id
                edgeid id = omega++;
                outEdgeIds[u][i] = id;
            }
        }
    });

    // copy edge ids for the edges in the other direction. Note that
    // "indexInOutEdgeArray" is slow which is why this second loop in parallel
    // makes sense.
    if (!directed) {
        balancedParallelForNodes([&](node u) {
            for (index i = 0; i < outEdges[u].size(); ++i) {
                node v = outEdges[u][i];
                if (v != none && outEdgeIds[u][i] == none) {
                    index j = indexInOutEdgeArray(v, u);
                    outEdgeIds[u][i] = outEdgeIds[v][j];
                }
            }
        });
    } else {
        balancedParallelForNodes([&](node u) {
            for (index i = 0; i < inEdges[u].size(); ++i) {
                node v = inEdges[u][i];
                if (v != none) {
                    index j = indexInOutEdgeArray(v, u);
                    inEdgeIds[u][i] = outEdgeIds[v][j];
                }
            }
        });
    }

    edgesIndexed = true; // remember that edges have been indexed so that addEdge
                         // needs to create edge ids
}

/** GRAPH INFORMATION **/

void GraphW::shrinkToFit() {
    exists.shrink_to_fit();

    inEdgeWeights.shrink_to_fit();
    for (auto &w : inEdgeWeights) {
        w.shrink_to_fit();
    }

    outEdgeWeights.shrink_to_fit();
    for (auto &w : outEdgeWeights) {
        w.shrink_to_fit();
    }

    inEdges.shrink_to_fit();
    for (auto &a : inEdges) {
        a.shrink_to_fit();
    }

    outEdges.shrink_to_fit();
    for (auto &a : outEdges) {
        a.shrink_to_fit();
    }
}

void GraphW::compactEdges() {
    this->parallelForNodes([&](node u) {
        if (degreeOut(u) == 0) {
            outEdges[u].clear();
            if (weighted)
                outEdgeWeights[u].clear();
            if (edgesIndexed)
                outEdgeIds[u].clear();
        } else {
            for (index i = 0; i < outEdges[u].size(); ++i) {
                while (i < outEdges[u].size() && outEdges[u][i] == none) {
                    outEdges[u][i] = outEdges[u].back();
                    outEdges[u].pop_back();

                    if (weighted) {
                        outEdgeWeights[u][i] = outEdgeWeights[u].back();
                        outEdgeWeights[u].pop_back();
                    }

                    if (edgesIndexed) {
                        outEdgeIds[u][i] = outEdgeIds[u].back();
                        outEdgeIds[u].pop_back();
                    }
                }
            }
        }
        if (directed) {
            if (degreeIn(u) == 0) {
                inEdges[u].clear();
                if (weighted)
                    inEdgeWeights[u].clear();
                if (edgesIndexed)
                    inEdgeIds[u].clear();
            } else {
                for (index i = 0; i < inEdges[u].size(); ++i) {
                    while (i < inEdges[u].size() && inEdges[u][i] == none) {
                        inEdges[u][i] = inEdges[u].back();
                        inEdges[u].pop_back();

                        if (weighted) {
                            inEdgeWeights[u][i] = inEdgeWeights[u].back();
                            inEdgeWeights[u].pop_back();
                        }

                        if (edgesIndexed) {
                            inEdgeIds[u][i] = inEdgeIds[u].back();
                            inEdgeIds[u].pop_back();
                        }
                    }
                }
            }
        }
    });
}

void GraphW::sortEdges() {
    std::vector<std::vector<node>> targetAdjacencies(upperNodeIdBound());
    std::vector<std::vector<edgeweight>> targetWeight;
    std::vector<std::vector<edgeid>> targetEdgeIds;

    if (isWeighted()) {
        targetWeight.resize(upperNodeIdBound());
        forNodes([&](node u) { targetWeight[u].reserve(degree(u)); });
    }
    if (hasEdgeIds()) {
        targetEdgeIds.resize(upperNodeIdBound());
        forNodes([&](node u) { targetEdgeIds[u].reserve(degree(u)); });
    }

    forNodes([&](node u) { targetAdjacencies[u].reserve(degree(u)); });

    auto assignToTarget = [&](node u, node v, edgeweight w, edgeid eid) {
        targetAdjacencies[v].push_back(u);
        if (isWeighted()) {
            targetWeight[v].push_back(w);
        }
        if (hasEdgeIds()) {
            targetEdgeIds[v].push_back(eid);
        }
    };

    forNodes([&](node u) { forInEdgesOf(u, assignToTarget); });

    outEdges.swap(targetAdjacencies);
    outEdgeWeights.swap(targetWeight);
    outEdgeIds.swap(targetEdgeIds);

    if (isDirected()) {
        inEdges.swap(targetAdjacencies);
        inEdgeWeights.swap(targetWeight);
        inEdgeIds.swap(targetEdgeIds);

        forNodes([&](node u) {
            targetAdjacencies[u].resize(degreeIn(u));
            targetAdjacencies[u].shrink_to_fit();
            targetAdjacencies[u].clear();
            if (isWeighted()) {
                targetWeight[u].resize(degreeIn(u));
                targetWeight[u].shrink_to_fit();
                targetWeight[u].clear();
            }
            if (hasEdgeIds()) {
                targetEdgeIds[u].resize(degreeIn(u));
                targetEdgeIds[u].shrink_to_fit();
                targetEdgeIds[u].clear();
            }
        });

        forNodes([&](node u) { forEdgesOf(u, assignToTarget); });

        inEdges.swap(targetAdjacencies);
        inEdgeWeights.swap(targetWeight);
        inEdgeIds.swap(targetEdgeIds);
    }
}

/** NODE MODIFIERS **/

node GraphW::addNode() {
    node v = z; // node gets maximum id
    z++;        // increment node range
    n++;        // increment node count

    // update per node data structures
    exists.push_back(true);

    outEdges.emplace_back();
    if (weighted)
        outEdgeWeights.emplace_back();
    if (edgesIndexed)
        outEdgeIds.emplace_back();

    if (directed) {
        inEdges.emplace_back();
        if (weighted)
            inEdgeWeights.emplace_back();
        if (edgesIndexed)
            inEdgeIds.emplace_back();
    }

    return v;
}

node GraphW::addNodes(count numberOfNewNodes) {
    if (numberOfNewNodes < 10) {
        // benchmarks suggested, it's cheaper to call 10 time emplace_back than resizing.
        while (numberOfNewNodes--)
            addNode();

        return z - 1;
    }

    z += numberOfNewNodes;
    n += numberOfNewNodes;

    // update per node data structures
    exists.resize(z, true);

    outEdges.resize(z);
    if (weighted)
        outEdgeWeights.resize(z);
    if (edgesIndexed)
        outEdgeIds.resize(z);

    if (directed) {
        inEdges.resize(z);
        if (weighted)
            inEdgeWeights.resize(z);
        if (edgesIndexed)
            inEdgeIds.resize(z);
    }

    return z - 1;
}

void GraphW::removeNode(node v) {
    assert(v < z);
    assert(exists[v]);

    // Remove all outgoing and ingoing edges
    while (!outEdges[v].empty())
        removeEdge(v, outEdges[v].front());
    if (isDirected())
        while (!inEdges[v].empty())
            removeEdge(inEdges[v].front(), v);

    // Make the attributes of this node invalid
    auto &theMap = nodeAttributes().attrMap;
    for (auto it = theMap.begin(); it != theMap.end(); ++it) {
        auto attributeStorageBase = it->second.get();
        attributeStorageBase->invalidate(v);
    }

    exists[v] = false;
    n--;
}

void GraphW::restoreNode(node v) {
    assert(v < z);
    assert(!exists[v]);

    exists[v] = true;
    n++;
}

/** EDGE MODIFIERS **/

bool GraphW::addEdge(node u, node v, edgeweight ew, bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && hasEdge(u, v)) {
        return false;
    }

    // increase number of edges
    ++m;
    outEdges[u].push_back(v);

    // if edges indexed, give new id
    if (edgesIndexed) {
        edgeid id = omega++;
        outEdgeIds[u].push_back(id);
    }

    if (directed) {
        inEdges[v].push_back(u);

        if (edgesIndexed) {
            inEdgeIds[v].push_back(omega - 1);
        }

        if (weighted) {
            inEdgeWeights[v].push_back(ew);
            outEdgeWeights[u].push_back(ew);
        }

    } else if (u == v) { // self-loop case
        if (weighted) {
            outEdgeWeights[u].push_back(ew);
        }
    } else { // undirected, no self-loop
        outEdges[v].push_back(u);

        if (weighted) {
            outEdgeWeights[u].push_back(ew);
            outEdgeWeights[v].push_back(ew);
        }

        if (edgesIndexed) {
            outEdgeIds[v].push_back(omega - 1);
        }
    }

    if (u == v) { // count self loop
        ++storedNumberOfSelfLoops;
    }

    return true;
}

bool GraphW::addPartialEdge(Unsafe, node u, node v, edgeweight ew, uint64_t index,
                           bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && (std::ranges::find(outEdges[u], v) != outEdges[u].end())) {
        return false;
    }

    outEdges[u].push_back(v);

    // if edges indexed, give new id
    if (edgesIndexed) {
        outEdgeIds[u].push_back(index);
    }
    if (weighted) {
        outEdgeWeights[u].push_back(ew);
    }

    return true;
}

bool GraphW::addPartialOutEdge(Unsafe, node u, node v, edgeweight ew, uint64_t index,
                              bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && (std::ranges::find(outEdges[u], v) != outEdges[u].end())) {
        return false;
    }

    outEdges[u].push_back(v);

    // if edges indexed, give new id
    if (edgesIndexed) {
        outEdgeIds[u].push_back(index);
    }
    if (weighted) {
        outEdgeWeights[u].push_back(ew);
    }

    return true;
}

bool GraphW::addPartialInEdge(Unsafe, node u, node v, edgeweight ew, uint64_t index,
                             bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && (std::ranges::find(inEdges[u], v) != inEdges[u].end())) {
        return false;
    }

    inEdges[u].push_back(v);

    if (edgesIndexed) {
        inEdgeIds[u].push_back(index);
    }
    if (weighted) {
        inEdgeWeights[u].push_back(ew);
    }

    return true;
}

template <typename T>
void erase(node u, index idx, std::vector<std::vector<T>> &vec) {
    vec[u][idx] = vec[u].back();
    vec[u].pop_back();
}

void GraphW::removeEdge(node u, node v) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (maintainCompactEdges && !edgesIndexed) {
        throw std::runtime_error("Edges have to be indexed if maintainCompactEdges is set to true");
    }

    index vi = indexInOutEdgeArray(u, v); // index in outEdges array
    index ui = indexInInEdgeArray(v, u);  // index in inEdges array

    if (edgesIndexed) {
        deletedID = edgeId(u, v);
    }

    if (vi == none) {
        std::stringstream strm;
        strm << "edge (" << u << "," << v << ") does not exist";
        throw std::runtime_error(strm.str());
    }

    const auto isLoop = (u == v);
    m--; // decrease number of edges
    if (isLoop)
        storedNumberOfSelfLoops--;

    // remove edge for source node
    erase<node>(u, vi, outEdges);
    if (weighted) {
        erase<edgeweight>(u, vi, outEdgeWeights);
    }
    if (edgesIndexed) {
        erase<edgeid>(u, vi, outEdgeIds);
        // Make the attributes of this edge invalid
        auto &theMap = edgeAttributes().attrMap;
        for (auto it = theMap.begin(); it != theMap.end(); ++it) {
            auto attributeStorageBase = it->second.get();
            attributeStorageBase->invalidate(deletedID);
        }
    }
    if (!directed && !isLoop) {
        // also remove edge for target node
        erase<node>(v, ui, outEdges);
        if (weighted) {
            erase<edgeweight>(v, ui, outEdgeWeights);
        }
        if (edgesIndexed) {
            erase<edgeid>(v, ui, outEdgeIds);
        }
    }
    if (maintainSortedEdges) {
        // initial index of deleted edge, also represents current index
        index cur = vi;

        // sort edges of source node from deleted index upwards
        while (cur + 1 < outEdges[u].size() && outEdges[u][cur] > outEdges[u][cur + 1]) {
            std::swap(outEdges[u][cur], outEdges[u][cur + 1]);
            if (edgesIndexed) {
                std::swap(outEdgeIds[u][cur], outEdgeIds[u][cur + 1]);
                // swap attributes as well
                auto &theMap = edgeAttributes().attrMap;
                for (auto it = theMap.begin(); it != theMap.end(); ++it) {
                    auto attributeStorageBase = it->second.get();
                    attributeStorageBase->swapData(outEdgeIds[u][cur], outEdgeIds[u][cur + 1]);
                }
            }
            ++cur;
        }

        if (!directed) {
            cur = ui;

            // sort edges of target node from deleted index upwards
            while (cur + 1 < outEdges[v].size() && outEdges[v][cur] > outEdges[v][cur + 1]) {
                std::swap(outEdges[v][cur], outEdges[v][cur + 1]);
                if (edgesIndexed) {
                    std::swap(outEdgeIds[v][cur], outEdgeIds[v][cur + 1]);
                }
                ++cur;
            }
        }
    }
    if (maintainCompactEdges) {
        // re-index edge IDs from deleted edge upwards
        balancedParallelForNodes([&](node w) {
            for (index i = 0; i < outEdges[w].size(); ++i) {
                auto curID = outEdgeIds[w][i];
                if (curID > deletedID) {
                    outEdgeIds[w][i]--;
                }
            }
        });
        // use erase to remove data entry at index `deletedID` and compact the data vector again
        auto &theMap = edgeAttributes().attrMap;
        for (auto it = theMap.begin(); it != theMap.end(); ++it) {
            auto attributeStorageBase = it->second.get();
            attributeStorageBase->erase(deletedID);
        }
    }
    if (directed) {
        assert(ui != none);

        erase<node>(v, ui, inEdges);
        if (weighted) {
            erase<edgeweight>(v, ui, inEdgeWeights);
        }
        if (edgesIndexed) {
            erase<edgeid>(v, ui, inEdgeIds);
        }
        if (maintainSortedEdges) {
            // initial index of deleted edge, also represents current index
            index cur = ui;

            // sort edges of target node from deleted index upwards
            while (cur + 1 < inEdges[v].size() && inEdges[v][cur] > inEdges[v][cur + 1]) {
                std::swap(inEdges[v][cur], inEdges[v][cur + 1]);
                if (edgesIndexed) {
                    std::swap(inEdgeIds[v][cur], inEdgeIds[v][cur + 1]);
                }
                ++cur;
            }
        }

        if (maintainCompactEdges) {
            // re-index edge ids from target node
            balancedParallelForNodes([&](node w) {
                for (index i = 0; i < inEdges[w].size(); ++i) {
                    node vv = inEdges[w][i];
                    if (vv != none) {
                        index j = indexInOutEdgeArray(vv, w);
                        inEdgeIds[w][i] = outEdgeIds[vv][j];
                    }
                }
            });
        }
    }
    if (maintainCompactEdges) {
        omega--; // decrease upperBound of edges
    }
}

void GraphW::removeAllEdges() {
    parallelForNodes([&](const node u) {
        removePartialOutEdges(unsafe, u);
        if (isDirected()) {
            removePartialInEdges(unsafe, u);
        }
    });

    m = 0;
}

void GraphW::removeSelfLoops() {
    parallelForNodes([&](const node u) {
        auto isSelfLoop = [u](const node v) { return u == v; };
        removeAdjacentEdges(u, isSelfLoop);
        if (isDirected()) {
            removeAdjacentEdges(u, isSelfLoop, true);
        }
    });

    m -= storedNumberOfSelfLoops;
    storedNumberOfSelfLoops = 0;
}

void GraphW::removeMultiEdges() {
    count removedEdges = 0;
    count removedSelfLoops = 0;
    std::unordered_set<node> nodes;

    forNodes([&](const node u) {
        nodes.reserve(degree(u));
        auto isMultiedge = [&nodes](const node v) { return !nodes.insert(v).second; };
        auto result = removeAdjacentEdges(u, isMultiedge);
        removedEdges += result.first;
        removedSelfLoops += result.second;
        if (isDirected()) {
            nodes.clear();
            removeAdjacentEdges(u, isMultiedge, true);
        }
        nodes.clear();
    });

    if (!isDirected()) {
        assert(!(removedEdges % 2));
        removedEdges /= 2;
    }

    m -= removedEdges + removedSelfLoops;
    storedNumberOfSelfLoops -= removedSelfLoops;
}

void GraphW::swapEdge(node s1, node t1, node s2, node t2) {
    index s1t1 = indexInOutEdgeArray(s1, t1);
    if (s1t1 == none)
        throw std::runtime_error("The first edge does not exist");
    index t1s1 = indexInInEdgeArray(t1, s1);

    index s2t2 = indexInOutEdgeArray(s2, t2);
    if (s2t2 == none)
        throw std::runtime_error("The second edge does not exist");
    index t2s2 = indexInInEdgeArray(t2, s2);

    std::swap(outEdges[s1][s1t1], outEdges[s2][s2t2]);

    if (directed) {
        std::swap(inEdges[t1][t1s1], inEdges[t2][t2s2]);

        if (weighted) {
            std::swap(inEdgeWeights[t1][t1s1], inEdgeWeights[t2][t2s2]);
        }

        if (edgesIndexed) {
            std::swap(inEdgeIds[t1][t1s1], inEdgeIds[t2][t2s2]);
        }
    } else {
        std::swap(outEdges[t1][t1s1], outEdges[t2][t2s2]);

        if (weighted) {
            std::swap(outEdgeWeights[t1][t1s1], outEdgeWeights[t2][t2s2]);
        }

        if (edgesIndexed) {
            std::swap(outEdgeIds[t1][t1s1], outEdgeIds[t2][t2s2]);
        }
    }
}

/** EDGE ATTRIBUTES **/

void GraphW::setWeight(node u, node v, edgeweight ew) {
    if (!weighted) {
        throw std::runtime_error("Cannot set edge weight in unweighted graph.");
    }

    index vi = indexInOutEdgeArray(u, v);
    if (vi == none) {
        // edge does not exist, create it, but warn user
        TRACE("Setting edge weight of a nonexisting edge will create the edge.");
        addEdge(u, v, ew);
        return;
    }

    outEdgeWeights[u][vi] = ew;
    if (directed) {
        index ui = indexInInEdgeArray(v, u);
        inEdgeWeights[v][ui] = ew;
    } else if (u != v) {
        index ui = indexInInEdgeArray(v, u);
        outEdgeWeights[v][ui] = ew;
    }
}

void GraphW::increaseWeight(node u, node v, edgeweight ew) {
    if (!weighted) {
        throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
    }

    index vi = indexInOutEdgeArray(u, v);
    if (vi == none) {
        // edge does not exits, create it, but warn user
        addEdge(u, v, ew);
        return;
    }

    outEdgeWeights[u][vi] += ew;
    if (directed) {
        index ui = indexInInEdgeArray(v, u);
        inEdgeWeights[v][ui] += ew;
    } else if (u != v) {
        index ui = indexInInEdgeArray(v, u);
        outEdgeWeights[v][ui] += ew;
    }
}

/** PREALLOCATION METHODS **/

void GraphW::preallocateUndirected(node u, size_t size) {
    assert(!directed);
    assert(exists[u]);
    outEdges[u].reserve(size);
    if (weighted) {
        outEdgeWeights[u].reserve(size);
    }
    if (edgesIndexed) {
        outEdgeIds[u].reserve(size);
    }
}

void GraphW::preallocateDirected(node u, size_t outSize, size_t inSize) {
    preallocateDirectedOutEdges(u, outSize);
    preallocateDirectedInEdges(u, inSize);
}

void GraphW::preallocateDirectedOutEdges(node u, size_t outSize) {
    assert(directed);
    assert(exists[u]);
    outEdges[u].reserve(outSize);

    if (weighted) {
        outEdgeWeights[u].reserve(outSize);
    }
    if (edgesIndexed) {
        outEdgeIds[u].reserve(outSize);
    }
}

void GraphW::preallocateDirectedInEdges(node u, size_t inSize) {
    assert(directed);
    assert(exists[u]);
    inEdges[u].reserve(inSize);

    if (weighted) {
        inEdgeWeights[u].reserve(inSize);
    }
    if (edgesIndexed) {
        inEdgeIds[u].reserve(inSize);
    }
}

} /* namespace NetworKit */
