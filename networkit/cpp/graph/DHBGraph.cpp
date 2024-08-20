#include <networkit/graph/DHBGraph.hpp>

#include <dhb/dynamic_hashed_blocks.h>

#include <cmath>
#include <map>
#include <random>
#include <sstream>
#include <thread>

namespace NetworKit {

/** CONSTRUCTORS **/

DHBGraph::DHBGraph(count n, bool weighted, bool directed, bool edgesIndexed)
    : n(n), m(0), storedNumberOfSelfLoops(0), z(n), omega(0),

      weighted(weighted),         // indicates whether the graph is weighted or not
      directed(directed),         // indicates whether the graph is directed or not
      edgesIndexed(edgesIndexed), // edges are not indexed by default

      exists(n, true),

      /* for directed graphs inEdges stores an adjacency list only considering
         incoming edges, for undirected graphs inEdges is not used*/
      inEdges(directed ? n : 0),

      /* for directed graphs outEdges stores an adjacency list only considering
      outgoing edges, for undirected graphs outEdges stores the adjacency list of
      undirected edges*/
      outEdges(n), inEdgeWeights(weighted && directed ? n : 0), outEdgeWeights(weighted ? n : 0),
      inEdgeIds(edgesIndexed && directed ? n : 0), outEdgeIds(edgesIndexed ? n : 0),
      nodeAttributeMap(this), edgeAttributeMap(this) {

#if defined(USING_DHB)
    m_dhb_graph = dhb::Matrix<EdgeData>(n);
#endif
}

DHBGraph::DHBGraph(std::initializer_list<WeightedEdge> edges) : DHBGraph(0, true) {
    using namespace std;

#if defined(USING_DHB)
    m_dhb_graph = dhb::Matrix<EdgeData>();
    for (const auto &weighted_edge : edges) {
        m_dhb_graph.insert(weighted_edge.u, weighted_edge.v, EdgeData{weighted_edge.weight, 0});
    }
#else

    /* Number of nodes = highest node index + 1 */
    for (const auto &edge : edges) {
        node x = std::max(edge.u, edge.v);
        while (numberOfNodes() <= x) {
            addNode();
        }
    }

    /* Now add all of the edges */
    for (const auto &edge : edges) {
        addEdge(edge.u, edge.v, edge.weight);
    }
#endif
}
/** PRIVATE HELPERS **/

index DHBGraph::indexInInEdgeArray(node v, node u) const {
    if (!directed) {
        return indexInOutEdgeArray(v, u);
    }
    for (index i = 0; i < inEdges[v].size(); i++) {
        node x = inEdges[v][i];
        if (x == u) {
            return i;
        }
    }
    return none;
}

index DHBGraph::indexInOutEdgeArray(node u, node v) const {
    for (index i = 0; i < outEdges[u].size(); i++) {
        node x = outEdges[u][i];
        if (x == v) {
            return i;
        }
    }
    return none;
}

/** EDGE IDS **/

void DHBGraph::indexEdges(bool force) {
    if (edgesIndexed && !force) {
        return;
    }
    omega = 0; // reset edge ids (for re-indexing)

#if defined(USING_DHB)
    if (directed) {
        forEdges([&](node u, node v) {
            m_dhb_graph.neighbors(u).iterator_to(v)->data().id = omega;
            ++omega;
        });
    } else {
        forEdges([&](node u, node v) {
            if (u == v) { // first time encounter self-loop
                m_dhb_graph.neighbors(u).iterator_to(v)->data().id = omega;
                ++omega;
            } else {
                m_dhb_graph.neighbors(u).iterator_to(v)->data().id = omega;
                m_dhb_graph.neighbors(v).iterator_to(u)->data().id = omega;
                ++omega;
            }
        });
    }
    edgesIndexed = true;

#else
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
#endif
}

edgeid DHBGraph::edgeId(node u, node v) const {
#if defined(USING_DHB)
    if (!edgesIndexed) {
        throw std::runtime_error(
            "Edges have not been indexed! To fix this, simply call indexEdges first.");
    }
    return m_dhb_graph.neighbors(u).iterator_to(v)->data().id;
#else

    if (!edgesIndexed) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    index i = indexInOutEdgeArray(u, v);

    if (i == none) {
        throw std::runtime_error("Edge does not exist");
    }
    edgeid id = outEdgeIds[u][i];
    return id;
#endif
}

/** GRAPH INFORMATION **/

bool DHBGraph::isIsolated(node v) const {
#if defined(USING_DHB)
    if (v >= m_dhb_graph.vertices_count()) {
        return true;
    }
    if (directed) {
        bool const has_no_outgoing_edges = (m_dhb_graph.degree(v) == 0);
        bool has_no_incoming_edges = true;
        forEdges([&](node from, node to) {
            if (v == to) {
                has_no_incoming_edges = false;
            }
        });
        return has_no_outgoing_edges && has_no_incoming_edges;
    }
    return (m_dhb_graph.degree(v) == 0);
#else
    if (!exists[v])
        throw std::runtime_error("Error, the node does not exist!");
    return outEdges[v].empty() && (!directed || inEdges[v].empty());
#endif
}

edgeweight DHBGraph::computeWeightedDegree(node u, bool inDegree, bool countSelfLoopsTwice) const {
    if (weighted) {
        edgeweight sum = 0.0;
        auto sumWeights = [&](node v, edgeweight w) {
            sum += (countSelfLoopsTwice && u == v) ? 2. * w : w;
        };
        if (inDegree) {
            forInNeighborsOf(u, sumWeights);
        } else {
            forNeighborsOf(u, sumWeights);
        }
        return sum;
    }

    count sum = inDegree ? degreeIn(u) : degreeOut(u);
    auto countSelfLoops = [&](node v) { sum += (u == v); };

    if (countSelfLoopsTwice && numberOfSelfLoops()) {
        if (inDegree) {
            forInNeighborsOf(u, countSelfLoops);
        } else {
            forNeighborsOf(u, countSelfLoops);
        }
    }

    return static_cast<edgeweight>(sum);
}

/** NODE MODIFIERS **/

node DHBGraph::addNode() {
#if defined(USING_DHB)
    m_dhb_graph.resize(m_dhb_graph.vertices_count() + 1);
    return m_dhb_graph.vertices_count() - 1;
#else
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
#endif
}

node DHBGraph::addNodes(count numberOfNewNodes) {
#if defined(USING_DHB)
    m_dhb_graph.resize(m_dhb_graph.vertices_count() + numberOfNewNodes);
    return m_dhb_graph.vertices_count() - 1;
#else
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
#endif
}

void DHBGraph::removeNode(node v) {
#if defined(USING_DHB)
    assert(v < m_dhb_graph.vertices_count());
    if (directed) {
        std::vector<std::pair<node, node>> edgesToRemove;
        forEdges([&](node from, node to) {
            if (from == v || to == v) {
                edgesToRemove.emplace_back(from, to);
            }
        });

        for (auto const& edge : edgesToRemove) {
            removeEdge(edge.first, edge.second);
        }
    } else {
        std::vector<node> neighbors;
        forNeighborsOf(v, [&](node u) { neighbors.push_back(u); });
        for (auto const& neighbor : neighbors) {
            removeEdge(v, neighbor);
            removeEdge(neighbor, v);
        }
    }

#else
    assert(v < z);
    assert(exists[v]);

    // Remove all outgoing and ingoing edges
    while (!outEdges[v].empty())
        removeEdge(v, outEdges[v].front());
    if (isDirected())
        while (!inEdges[v].empty())
            removeEdge(inEdges[v].front(), v);

    // Make the attributes of this node invalid
    auto& theMap = nodeAttributeMap.attrMap;
    for (auto it = theMap.begin(); it != theMap.end(); ++it) {
        auto attributeStorageBase = it->second.get();
        attributeStorageBase->invalidate(v);
    }

    exists[v] = false;
    n--;
#endif
}

/** NODE PROPERTIES **/

edgeweight DHBGraph::weightedDegree(node u, bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, false, countSelfLoopsTwice);
}

edgeweight DHBGraph::weightedDegreeIn(node u, bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, true, countSelfLoopsTwice);
}

bool DHBGraph::hasNode(node v) const noexcept {
#if defined(USING_DHB)
    return (v < m_dhb_graph.vertices_count());
#else
    return (v < z) && this->exists[v];
#endif
}

/** EDGE MODIFIERS **/

bool DHBGraph::addEdge(node u, node v, edgeweight ew, bool checkForMultiEdges) {
#if defined(USING_DHB)
    assert(u < m_dhb_graph.vertices_count());
    assert(v < m_dhb_graph.vertices_count());

    edgeid const id = edgesIndexed ? omega : 0;

    auto const edge_weight = weighted ? ew : defaultEdgeWeight;
    EdgeData const edge_data{edge_weight, dhb::EdgeID{id}};

    bool is_inserted_in = false;
    bool is_inserted_out = false;
    bool is_inserted = false;

    // Insert edges depending on whether the graph is directed or undirected
    if (!directed) {
        is_inserted_out = m_dhb_graph.insert(u, v, edge_data);
        is_inserted_in = m_dhb_graph.insert(v, u, edge_data);
    } else {
        is_inserted = m_dhb_graph.insert(u, v, edge_data);
    }

    bool const undirected_insertion_success = is_inserted_in && is_inserted_out && !directed;
    bool const directed_insertion_success = is_inserted && directed;

    if (undirected_insertion_success || directed_insertion_success) {
        if (u == v) {
            ++storedNumberOfSelfLoops;
        }
        if (edgesIndexed) {
            omega++;
        }
        return true;
    }
    return false;
#else
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
#endif
}

bool DHBGraph::addEdges(std::vector<WeightedEdge>&& weighted_edges, bool do_update,
                        unsigned int num_threads) {
#if defined(USING_DHB)
    omp_set_num_threads(num_threads);
    assert(omp_get_max_threads() == num_threads);

    auto cmp = [](WeightedEdge const& a, WeightedEdge const& b) { return a.u < b.u; };
    auto key = [](WeightedEdge const& e) { return e.u; };

    std::vector<uint8_t> insertion_result(omp_get_max_threads(), 1);
    auto insert_edge_f = [&](WeightedEdge const& e) {
        EdgeData data{e.weight, 0};
        auto [it, inserted] = m_dhb_graph.neighbors(e.u).insert(e.v, data);
        if (inserted) {
            storedNumberOfSelfLoops += uint64_t(e.u == e.v);
        }
        uint8_t& insertion_result_t = insertion_result[omp_get_thread_num()];
        insertion_result_t = insertion_result_t && inserted;
    };

    auto insert_and_update_f = [&](WeightedEdge const& e) {
        EdgeData data{e.weight, 0};
        auto [it, inserted] = m_dhb_graph.neighbors(e.u).insert(e.v, data);
        if (inserted) {
            storedNumberOfSelfLoops += uint64_t(e.u == e.v);
        } else {
            it->data() = EdgeData{e.weight, it->data().id};
        }
        uint8_t& insertion_result_t = insertion_result[omp_get_thread_num()];
        insertion_result_t = insertion_result_t && inserted;
    };

    auto get_source_f = [](WeightedEdge const& e) { return e.u; };
    dhb::BatchParallelizer<WeightedEdge> par;
    auto processEdge = [&](WeightedEdge const& e) {
        if (do_update) {
            insert_and_update_f(e);
        } else {
            insert_edge_f(e);
        }
    };

    par(std::begin(weighted_edges), std::end(weighted_edges), std::move(get_source_f),
        std::move(key), std::move(cmp), std::move(processEdge));

    bool const acc_insertion_result_directed = std::all_of(
        std::begin(insertion_result), std::end(insertion_result), [](bool const r) { return r; });

    bool const acc_insertion_result_undirected = [&]() {
        if (!isDirected()) {
            std::vector<WeightedEdge> weighted_edges_v_to_u;
            weighted_edges_v_to_u.resize(weighted_edges.size());

#pragma omp parallel for
            for (count i = 0; i < weighted_edges_v_to_u.size(); ++i) {
                weighted_edges_v_to_u[i].u = weighted_edges[i].v;
                weighted_edges_v_to_u[i].v = weighted_edges[i].u;
                weighted_edges_v_to_u[i].weight = weighted_edges[i].weight;
            }

            par(std::begin(weighted_edges_v_to_u), std::end(weighted_edges_v_to_u),
                std::move(get_source_f), std::move(key), std::move(cmp), std::move(processEdge));

            return std::all_of(std::begin(insertion_result), std::end(insertion_result),
                               [](bool const r) { return r; });
        }

        return true;
    }();

    return acc_insertion_result_directed && acc_insertion_result_undirected;
#endif
    return false;
}

bool DHBGraph::addEdges(std::vector<Edge>&& edges, bool do_update, unsigned int num_threads) {
    std::vector<WeightedEdge> weighted_edges;
    weighted_edges.reserve(edges.size());
    for (auto& edge : edges) {
        WeightedEdge w_edge = WeightedEdge{edge.u, edge.v, defaultEdgeWeight};
        weighted_edges.emplace_back(std::move(w_edge));
    }
    return addEdges(std::move(weighted_edges), do_update, num_threads);
}

template <typename T>
void erase(node u, index idx, std::vector<std::vector<T>> &vec) {
    vec[u][idx] = vec[u].back();
    vec[u].pop_back();
}

bool DHBGraph::removeEdge(node u, node v) {
#if defined(USING_DHB)
    bool const is_loop = (u == v);
    bool is_removed_in = false;
    bool is_removed_out = false;
    bool is_removed = false;

    if (!directed) {
        is_removed_out = m_dhb_graph.removeEdge(u, v);
        is_removed_in = m_dhb_graph.removeEdge(v, u);
    } else {
        is_removed = m_dhb_graph.removeEdge(u, v);
    }

    bool const undirected_removal_success = (is_removed_in || is_removed_out) && !isDirected();
    bool const directed_removal_success = is_removed && isDirected();

    if (undirected_removal_success || directed_removal_success) {
        if (is_loop) {
            storedNumberOfSelfLoops--;
        }
        return true;
    }
    return false;
#else
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (maintainCompactEdges && !edgesIndexed) {
        throw std::runtime_error("Edges have to be indexed if maintainCompactEdges is set to true");
    }

    index vi = indexInOutEdgeArray(u, v);
    index ui = indexInInEdgeArray(v, u);

    if (maintainCompactEdges) {
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
        // Make the attributes of this node invalid
        auto &theMap = nodeAttributeMap.attrMap;
        for (auto it = theMap.begin(); it != theMap.end(); ++it) {
            auto attributeStorageBase = it->second.get();
            attributeStorageBase->invalidate(vi);
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
        balancedParallelForNodes([&](node u) {
            for (index i = 0; i < outEdges[u].size(); ++i) {
                auto curID = outEdgeIds[u][i];
                if (curID > deletedID) {
                    outEdgeIds[u][i]--;
                }
            }
        });
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
    }
    if (maintainCompactEdges) {
        omega--; // decrease upperBound of edges
    }
    return true;
#endif
}

void DHBGraph::removeAllEdges() {
#if defined(USING_DHB)
    m_dhb_graph = dhb::Matrix<EdgeData>(m_dhb_graph.vertices_count());
#else
    forNodesParallel([&](node const u) {
        removePartialOutEdges(unsafe, u);
        if (isDirected()) {
            removePartialInEdges(unsafe, u);
        }
    });

    m = 0;
#endif
}

void DHBGraph::removeSelfLoops() {
    forNodesParallel([&](node const u) {
        auto isSelfLoop = [u](node const v) { return u == v; };
        removeAdjacentEdges(u, isSelfLoop);
        if (isDirected()) {
            removeAdjacentEdges(u, isSelfLoop, true);
        }
    });

    m -= storedNumberOfSelfLoops;
    storedNumberOfSelfLoops = 0;
}

count DHBGraph::numberOfEdges() const noexcept {
#if defined(USING_DHB)
    // DHB stores only directed edges. Thus, if the count of edges for an undirected graph
    // is requested we must divide the count of edges coming from the DHB object by 2.

    // DHB stores the edge only once if a self loop is added, because it does not allow duplicate
    // edges. Therefore storedNumberOfSelfLoops should be added to count for counting edges
    // correctly.
    uint64_t count = m_dhb_graph.edges_count();
    if (directed) {
        return count;
    }
    count += storedNumberOfSelfLoops;
    assert(count % 2 == 0);
    return count >> 1;
#else
    return m;
#endif
}

index DHBGraph::upperNodeIdBound() const noexcept {
#if defined(USING_DHB)
    return m_dhb_graph.vertices_count();
#else
    return z;
#endif
}

count DHBGraph::numberOfNodes() const noexcept {
#if defined(USING_DHB)
    return m_dhb_graph.vertices_count();
#else
    return n;
#endif
}

count DHBGraph::degree(node v) const {
#if defined(USING_DHB)
    assert(v < m_dhb_graph.vertices_count());
    return m_dhb_graph.degree(v);
#else
    assert(hasNode(v));
    return outEdges[v].size();
#endif
}

count DHBGraph::degreeIn(node v) const {
#if defined(USING_DHB)
    assert(v < m_dhb_graph.vertices_count());
    if (directed) {
        count in_degree = 0;
        parallelForEdges([&](node from, node to) {
            if (v == to) {
#pragma omp atomic
                ++in_degree;
            }
        });
        return in_degree;
    }
    return degreeOut(v);

#else
    assert(hasNode(v));
    return directed ? inEdges[v].size() : outEdges[v].size();
#endif
}

void DHBGraph::swapEdge(node s1, node t1, node s2, node t2) {
#if defined(USING_DHB)
    edgeweight const weight_1 = isWeighted() ? weight(s1, t1) : defaultEdgeWeight;
    edgeweight const weight_2 = isWeighted() ? weight(s2, t2) : defaultEdgeWeight;
    edgeid const id_1 = edgeId(s1, t1);
    edgeid const id_2 = edgeId(s2, t2);

    if (isDirected()) {
        m_dhb_graph.neighbors(s1).iterator_to(t1).get_ptr()->vertex = t2;
        auto iter_data_t1 = m_dhb_graph.neighbors(s1).iterator_to(t1)->data();
        iter_data_t1.weight = weight_2;
        iter_data_t1.id = id_2;

        m_dhb_graph.neighbors(s2).iterator_to(t2).get_ptr()->vertex = t1;
        auto iter_data_t2 = m_dhb_graph.neighbors(s2).iterator_to(t2)->data();
        iter_data_t2.weight = weight_1;
        iter_data_t2.id = id_1;
    } else {
        removeEdge(s1, t1);
        removeEdge(s2, t2);
        m_dhb_graph.insert(s1, t2, EdgeData{weight_1, id_1});
        m_dhb_graph.insert(t2, s1, EdgeData{weight_1, id_1});
        m_dhb_graph.insert(s2, t1, EdgeData{weight_2, id_2});
        m_dhb_graph.insert(t1, s2, EdgeData{weight_2, id_2});
    }
    if (t1 == s2 || t2 == s1) {
        storedNumberOfSelfLoops++;
    }

#else
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
#endif
}

bool DHBGraph::hasEdge(node u, node v) const noexcept {
#if defined(USING_DHB)
    return m_dhb_graph.neighbors(u).exists(v);

#else
    if (u >= z || v >= z) {
        return false;
    }
    if (!directed && outEdges[u].size() > outEdges[v].size()) {
        return indexInOutEdgeArray(v, u) != none;
    } else if (directed && outEdges[u].size() > inEdges[v].size()) {
        return indexInInEdgeArray(v, u) != none;
    } else {
        return indexInOutEdgeArray(u, v) != none;
    }
#endif
}

/** EDGE ATTRIBUTES **/

edgeweight DHBGraph::weight(node u, node v) const {
#if defined(USING_DHB)
    auto it = m_dhb_graph.neighbors(u).iterator_to(v);
    if (it != m_dhb_graph.neighbors(u).end()) {
        return isWeighted() ? it->data().weight : defaultEdgeWeight;
    }
    return nullWeight;
#else

    index vi = indexInOutEdgeArray(u, v);
    if (vi == none) {
        return nullWeight;
    } else {
        return weighted ? outEdgeWeights[u][vi] : defaultEdgeWeight;
    }
#endif
}

void DHBGraph::setWeight(node u, node v, edgeweight ew) {
#if defined(USING_DHB)
    if (!hasEdge(u, v)) {
        // edge does not exist, create it, but warn user
        TRACE("Setting edge weight of a nonexisting edge will create the edge.");
        addEdge(u, v, ew);
        return;
    }

    m_dhb_graph.neighbors(u).iterator_to(v)->data().weight = isWeighted() ? ew : defaultEdgeWeight;
    if (!directed) {
        m_dhb_graph.neighbors(v).iterator_to(u)->data().weight =
            isWeighted() ? ew : defaultEdgeWeight;
    }
#else

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
#endif
}

void DHBGraph::increaseWeight(node u, node v, edgeweight ew) {
    if (!weighted) {
        throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
    }
#if defined(USING_DHB)
    if (!hasEdge(u, v)) {
        // edge does not exist, create it, but warn user
        TRACE("Increasing the weight of a nonexistent edge will create the edge with the proposed "
              "increased weight.");
        addEdge(u, v, ew);
        return;
    }
    edgeweight const new_ew = ew + weight(u, v);
    m_dhb_graph.neighbors(u).iterator_to(v)->data().weight = new_ew;
    if (!directed) {
        m_dhb_graph.neighbors(v).iterator_to(u)->data().weight = new_ew;
    }
#else
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
#endif
}

void DHBGraph::setWeightAtIthNeighbor(Unsafe, node u, index i, edgeweight ew) {
#if defined(USING_DHB)
    m_dhb_graph.neighbors(u)[i].data().weight = ew;
#else
    outEdgeWeights[u][i] = ew;
#endif
}

void DHBGraph::setWeightAtIthInNeighbor(Unsafe, node u, index i, edgeweight ew) {
#if defined(USING_DHB)
    setWeight(getIthInNeighbor(u, i), u, ew);
#else
    inEdgeWeights[u][i] = ew;
#endif
}

index DHBGraph::indexOfNeighbor(node u, node v) const {
#if defined(USING_DHB)
    for (index i = 0; i < m_dhb_graph.neighbors(u).degree(); i++) {
        node current = m_dhb_graph.neighbors(u)[i].vertex();
        if (v == current) {
            return i;
        }
    }
    return none;
#else
    return indexInOutEdgeArray(u, v);
#endif
}

node DHBGraph::getIthNeighbor(node u, index i) const {
#if defined(USING_DHB)
    if (!ithNeighborExists(u, i)) {
        return none;
    }
    return m_dhb_graph.neighbors(u)[i].vertex();
#else
    if (!hasNode(u) || i >= outEdges[u].size())
        return none;
    return outEdges[u][i];
#endif
}

node DHBGraph::getIthInNeighbor(node u, index i) const {
#if defined(USING_DHB)
    index currentIndex = 0;
    for (dhb::Vertex from = 0; from < m_dhb_graph.vertices_count(); ++from) {
        for (auto it = neighborRange(from).begin(); it != neighborRange(from).end(); ++it) {
            dhb::Vertex const to = *it;
            if (u == to) {
                if (currentIndex == i) {
                    return from;
                }
                currentIndex++;
            }
        }
    }
    return none;
#else
    if (!hasNode(u) || i >= inEdges[u].size())
        return none;
    return inEdges[u][i];
#endif
}

edgeweight DHBGraph::getIthNeighborWeight(node u, index i) const {
#if defined(USING_DHB)
    if (!ithNeighborExists(u, i)) {
        return none;
    }
    return isWeighted() ? m_dhb_graph.neighbors(u)[i].data().weight : defaultEdgeWeight;
#else
    if (!hasNode(u) || i >= outEdges[u].size())
        return nullWeight;
    return isWeighted() ? outEdgeWeights[u][i] : defaultEdgeWeight;
#endif
}

std::pair<node, edgeweight> DHBGraph::getIthNeighborWithWeight(node u, index i) const {
#if defined(USING_DHB)
    if (!ithNeighborExists(u, i)) {
        return {none, none};
    }
    node const neighbor = m_dhb_graph.neighbors(u)[i].vertex();
    double const weight =
        isWeighted() ? m_dhb_graph.neighbors(u)[i].data().weight : defaultEdgeWeight;
    return {neighbor, weight};
#else
    if (!hasNode(u) || i >= outEdges[u].size())
        return {none, none};
    return getIthNeighborWithWeight(unsafe, u, i);
#endif
}

std::pair<node, edgeweight> DHBGraph::getIthNeighborWithWeight(Unsafe, node u, index i) const {
#if defined(USING_DHB)
    return {none, none};
#else
    if (!isWeighted())
        return {outEdges[u][i], defaultEdgeWeight};
    return {outEdges[u][i], outEdgeWeights[u][i]};
#endif
}

std::pair<node, edgeid> DHBGraph::getIthNeighborWithId(node u, index i) const {
    assert(hasEdgeIds());
#if defined(USING_DHB)
    if (!ithNeighborExists(u, i)) {
        return {none, none};
    }
    node const neighbor = m_dhb_graph.neighbors(u)[i].vertex();
    edgeid const id = m_dhb_graph.neighbors(u)[i].data().id;
    return {neighbor, id};
#else
    if (!hasNode(u) || i >= outEdges[u].size())
        return {none, none};
    return {outEdges[u][i], outEdgeIds[u][i]};
#endif
}

/** SUMS **/

edgeweight DHBGraph::totalEdgeWeight() const noexcept {
    if (weighted)
        return sumForEdgesParallel([](node, node, edgeweight ew) { return ew; });
    else
        return numberOfEdges() * defaultEdgeWeight;
}

/* ATTRIBUTE PREMISE AND INDEX CHECKS */

template <>
void DHBGraph::ASB<DHBGraph::PerNode>::checkPremise() const {
    // nothing
}

template <>
void DHBGraph::ASB<DHBGraph::PerEdge>::checkPremise() const {
    if (!theGraph->hasEdgeIds()) {
        throw std::runtime_error("Edges must be indexed");
    }
}

template <>
void DHBGraph::ASB<DHBGraph::PerNode>::indexOK(index n) const {
    if (!theGraph->hasNode(n)) {
        throw std::runtime_error("This node does not exist");
    }
}

template <>
void DHBGraph::ASB<DHBGraph::PerEdge>::indexOK(index n) const {
    auto uv = theGraph->edgeById(n);
    if (!theGraph->hasEdge(uv.first, uv.second)) {
        throw std::runtime_error("This edgeId does not exist");
    }
}

} // namespace NetworKit
