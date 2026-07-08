#ifndef NETWORKIT_GRAPH_ADJ_LIST_GRAPH_IMPL_HPP_
#define NETWORKIT_GRAPH_ADJ_LIST_GRAPH_IMPL_HPP_

namespace NetworKit {

/* NODE ITERATORS */

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    for (NodeT v = 0; v < z; ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::parallelForNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
#pragma omp parallel for
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename C, typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodesWhile(
    C condition, L handle) const { // NOLINT(performance-unnecessary-value-param)
    for (NodeT v = 0; v < z; ++v) {
        if (exists[v]) {
            if (!condition()) {
                break;
            }
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodesInRandomOrder(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    std::vector<NodeT> randVec;
    randVec.reserve(numberOfNodes());
    forNodes([&](NodeT u) { randVec.push_back(u); });
    std::ranges::shuffle(randVec, Aux::Random::getURNG());
    for (NodeT v : randVec) {
        handle(v);
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::balancedParallelForNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
// TODO: define min block size (and test it!)
#pragma omp parallel for schedule(guided)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodePairs(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    for (NodeT u = 0; u < z; ++u) {
        if (exists[u]) {
            for (NodeT v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::parallelForNodePairs(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        if (exists[u]) {
            for (NodeT v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

/* EDGE ITERATORS */

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::forOutEdgesOfImpl(NodeT u, L handle) const {
    for (index i = 0; i < outEdges[u].size(); ++i) {
        NodeT v = outEdges[u][i];

        if (useEdgeInIteration<graphIsDirected>(u, v)) {
            edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                          getOutEdgeId<graphHasEdgeIds>(u, i));
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::forInEdgesOfImpl(NodeT u, L handle) const {
    if (graphIsDirected) {
        for (index i = 0; i < inEdges[u].size(); ++i) {
            NodeT v = inEdges[u][i];

            edgeLambda<L>(handle, u, v, getInEdgeWeight<hasWeights>(u, i),
                          getInEdgeId<graphHasEdgeIds>(u, i));
        }
    } else {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            NodeT v = outEdges[u][i];

            edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                          getOutEdgeId<graphHasEdgeIds>(u, i));
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::forEdgeImpl(L handle) const {
    for (NodeT u = 0; u < z; ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(u, handle);
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::parallelForEdgesImpl(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(u, handle);
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline double AdjListGraph<NodeT, EdgeWeightT>::parallelSumForEdgesImpl(L handle) const {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            NodeT v = outEdges[u][i];

            // undirected, do not iterate over edges twice
            // {u, v} instead of (u, v); if v == nullNodeId, u > v is not fulfilled
            if (useEdgeInIteration<graphIsDirected>(u, v)) {
                sum += edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                                     getOutEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    }

    return sum;
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forEdges(L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edgeIds
        forEdgeImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edgeIds
        forEdgeImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edgeIds
        forEdgeImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted, directed, no edgeIds
        forEdgeImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edgeIds
        forEdgeImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edgeIds
        forEdgeImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edgeIds
        forEdgeImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edgeIds
        forEdgeImpl<true, true, true, L>(handle);
        break;
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::parallelForEdges(L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edgeIds
        parallelForEdgesImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edgeIds
        parallelForEdgesImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edgeIds
        parallelForEdgesImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted, directed, no edgeIds
        parallelForEdgesImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edgeIds
        parallelForEdgesImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edgeIds
        parallelForEdgesImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edgeIds
        parallelForEdgesImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edgeIds
        parallelForEdgesImpl<true, true, true, L>(handle);
        break;
    }
}

/* NEIGHBORHOOD ITERATORS */

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNeighborsOf(NodeT u, L handle) const {
    forEdgesOf(u, handle);
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forEdgesOf(NodeT u, L handle) const {
    switch (weighted + 2 * edgesIndexed) {
    case 0: // not weighted, no edge ids
        forOutEdgesOfImpl<true, false, false, L>(u, handle);
        break;

    case 1: // weighted, no edge ids
        forOutEdgesOfImpl<true, true, false, L>(u, handle);
        break;

    case 2: // not weighted, with edge ids
        forOutEdgesOfImpl<true, false, true, L>(u, handle);
        break;

    case 3: // weighted, with edge ids
        forOutEdgesOfImpl<true, true, true, L>(u, handle);
        break;
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forInNeighborsOf(NodeT u, L handle) const {
    forInEdgesOf(u, handle);
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forInEdgesOf(NodeT u, L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edge ids
        forInEdgesOfImpl<false, false, false, L>(u, handle);
        break;

    case 1: // weighted, undirected, no edge ids
        forInEdgesOfImpl<false, true, false, L>(u, handle);
        break;

    case 2: // unweighted, directed, no edge ids
        forInEdgesOfImpl<true, false, false, L>(u, handle);
        break;

    case 3: // weighted, directed, no edge ids
        forInEdgesOfImpl<true, true, false, L>(u, handle);
        break;

    case 4: // unweighted, undirected, with edge ids
        forInEdgesOfImpl<false, false, true, L>(u, handle);
        break;

    case 5: // weighted, undirected, with edge ids
        forInEdgesOfImpl<false, true, true, L>(u, handle);
        break;

    case 6: // unweighted, directed, with edge ids
        forInEdgesOfImpl<true, false, true, L>(u, handle);
        break;

    case 7: // weighted, directed, with edge ids
        forInEdgesOfImpl<true, true, true, L>(u, handle);
        break;
    }
}

/* REDUCTION ITERATORS */

template <class NodeT, class EdgeWeightT>
template <typename L>
double AdjListGraph<NodeT, EdgeWeightT>::parallelSumForNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            sum += handle(v);
        }
    }

    return sum;
}

template <class NodeT, class EdgeWeightT>
template <typename L>
double AdjListGraph<NodeT, EdgeWeightT>::parallelSumForEdges(L handle) const {
    double sum = 0.0;

    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edge ids
        sum = parallelSumForEdgesImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edge ids
        sum = parallelSumForEdgesImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edge ids
        sum = parallelSumForEdgesImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted,   directed, no edge ids
        sum = parallelSumForEdgesImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edge ids
        sum = parallelSumForEdgesImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edge ids
        sum = parallelSumForEdgesImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edge ids
        sum = parallelSumForEdgesImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edge ids
        sum = parallelSumForEdgesImpl<true, true, true, L>(handle);
        break;
    }

    return sum;
}

/* EDGE MODIFIERS */

template <class NodeT, class EdgeWeightT>
template <typename Condition>
std::pair<count, count>
AdjListGraph<NodeT, EdgeWeightT>::removeAdjacentEdges(NodeT u, Condition condition, bool edgesIn) {
    count removedEdges = 0;
    count removedSelfLoops = 0;

    // For directed graphs, this function is supposed to be called twice: one to remove out-edges,
    // and one to remove in-edges.
    auto &edges_ = edgesIn ? inEdges[u] : outEdges[u];
    for (index vi = 0; vi < edges_.size();) {
        if (condition(edges_[vi])) {
            const auto isSelfLoop = (edges_[vi] == u);
            removedSelfLoops += isSelfLoop;
            removedEdges += !isSelfLoop;
            edges_[vi] = edges_.back();
            edges_.pop_back();
            if (isWeighted()) {
                auto &weights_ = edgesIn ? inEdgeWeights[u] : outEdgeWeights[u];
                weights_[vi] = weights_.back();
                weights_.pop_back();
            }
            if (hasEdgeIds()) {
                auto &edgeIds_ = edgesIn ? inEdgeIds[u] : outEdgeIds[u];
                edgeIds_[vi] = edgeIds_.back();
                edgeIds_.pop_back();
            }
        } else {
            ++vi;
        }
    }

    return {removedEdges, removedSelfLoops};
}

template <class NodeT, class EdgeWeightT>
template <typename Lambda>
void AdjListGraph<NodeT, EdgeWeightT>::sortNeighbors(NodeT u, Lambda lambda) {
    if ((degreeIn(u) < 2) && (degree(u) < 2)) {
        return;
    }
    // Sort the outEdge-Attributes
    std::vector<index> outIndices(outEdges[u].size());
    std::iota(outIndices.begin(), outIndices.end(), 0);
    std::ranges::sort(outIndices,
                      [&](index a, index b) { return lambda(outEdges[u][a], outEdges[u][b]); });

    Aux::ArrayTools::applyPermutation(outEdges[u].begin(), outEdges[u].end(), outIndices.begin());

    if (weighted) {
        Aux::ArrayTools::applyPermutation(outEdgeWeights[u].begin(), outEdgeWeights[u].end(),
                                          outIndices.begin());
    }

    if (edgesIndexed) {
        Aux::ArrayTools::applyPermutation(outEdgeIds[u].begin(), outEdgeIds[u].end(),
                                          outIndices.begin());
    }

    // For directed graphs we need to sort the inEdge-Attributes separately
    if (directed) {
        std::vector<index> inIndices(inEdges[u].size());
        std::iota(inIndices.begin(), inIndices.end(), 0);

        std::ranges::sort(inIndices,
                          [&](index a, index b) { return lambda(inEdges[u][a], inEdges[u][b]); });

        Aux::ArrayTools::applyPermutation(inEdges[u].begin(), inEdges[u].end(), inIndices.begin());

        if (weighted) {
            Aux::ArrayTools::applyPermutation(inEdgeWeights[u].begin(), inEdgeWeights[u].end(),
                                              inIndices.begin());
        }

        if (edgesIndexed) {
            Aux::ArrayTools::applyPermutation(inEdgeIds[u].begin(), inEdgeIds[u].end(),
                                              inIndices.begin());
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <class Lambda>
void AdjListGraph<NodeT, EdgeWeightT>::sortEdges(Lambda lambda) {

    std::vector<std::vector<index>> indicesGlobal(omp_get_max_threads());

    const auto sortAdjacencyArrays = [&](NodeT u, std::vector<NodeT> &adjList,
                                         std::vector<EdgeWeightT> &weights,
                                         std::vector<edgeid> &edgeIds) -> void {
        auto &indices = indicesGlobal[omp_get_thread_num()];
        if (adjList.size() > indices.size())
            indices.resize(adjList.size());

        const auto indicesEnd =
            indices.begin()
            + static_cast<
                std::iterator_traits<std::vector<index>::const_iterator>::difference_type>(
                adjList.size());
        std::iota(indices.begin(), indicesEnd, 0);

        if (isWeighted()) {
            if (hasEdgeIds())
                std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                    return lambda(WeightedEdgeWithId{u, adjList[a], weights[a], edgeIds[a]},
                                  WeightedEdgeWithId{u, adjList[b], weights[b], edgeIds[b]});
                });
            else
                std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                    return lambda(WeightedEdgeWithId{u, adjList[a], weights[a], 0},
                                  WeightedEdgeWithId{u, adjList[b], weights[b], 0});
                });
        } else if (hasEdgeIds())
            std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                return lambda(WeightedEdgeWithId{u, adjList[a], defaultEdgeWeightT, edgeIds[a]},
                              WeightedEdgeWithId{u, adjList[b], defaultEdgeWeightT, edgeIds[b]});
            });
        else
            std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                return lambda(WeightedEdgeWithId{u, adjList[a], defaultEdgeWeightT, 0},
                              WeightedEdgeWithId{u, adjList[b], defaultEdgeWeightT, 0});
            });

        Aux::ArrayTools::applyPermutation(adjList.begin(), adjList.end(), indices.begin());

        if (isWeighted())
            Aux::ArrayTools::applyPermutation(weights.begin(), weights.end(), indices.begin());

        if (hasEdgeIds())
            Aux::ArrayTools::applyPermutation(edgeIds.begin(), edgeIds.end(), indices.begin());
    };

    balancedParallelForNodes([&](const NodeT u) {
        if (degree(u) < 2)
            return;

        std::vector<EdgeWeightT> dummyEdgeWeights;
        std::vector<edgeid> dummyEdgeIds;
        sortAdjacencyArrays(u, outEdges[u], isWeighted() ? outEdgeWeights[u] : dummyEdgeWeights,
                            hasEdgeIds() ? outEdgeIds[u] : dummyEdgeIds);

        if (isDirected())
            sortAdjacencyArrays(u, inEdges[u], isWeighted() ? inEdgeWeights[u] : dummyEdgeWeights,
                                hasEdgeIds() ? inEdgeIds[u] : dummyEdgeIds);
    });
}

/** CONSTRUCTORS **/

template <class NodeT, class EdgeWeightT>
AdjListGraph<NodeT, EdgeWeightT>::AdjListGraph(count n, bool weighted, bool directed,
                                               bool edgesIndexed)
    : n(n), m(0), storedNumberOfSelfLoops(0), z(n), omega(0), t(0),

      weighted(weighted), // indicates whether the graph is weighted or not
      directed(directed), // indicates whether the graph is directed or not
      edgesIndexed(edgesIndexed), deletedID(nullEdgeId),
      // edges are not indexed by default

      exists(n, true),

      /* for directed graphs inEdges stores an adjacency list only considering
         incoming edges, for undirected graphs inEdges is not used*/
      inEdges(directed ? n : 0),

      /* for directed graphs outEdges stores an adjacency list only considering
      outgoing edges, for undirected graphs outEdges stores the adjacency list of
      undirected edges*/
      outEdges(n), inEdgeWeights(weighted && directed ? n : 0), outEdgeWeights(weighted ? n : 0),
      inEdgeIds(edgesIndexed && directed ? n : 0), outEdgeIds(edgesIndexed ? n : 0),
      nodeAttributeMap(this), edgeAttributeMap(this) {}

template <class NodeT, class EdgeWeightT>
AdjListGraph<NodeT, EdgeWeightT>::AdjListGraph(
    std::initializer_list<WeightedEdgeT<NodeT, EdgeWeightT>> edges)
    : AdjListGraph(0, true) {

    /* Number of nodes = highest node index + 1 */
    for (const auto &edge : edges) {
        NodeT x = std::max(edge.u, edge.v);
        while (numberOfNodes() <= x) {
            addNode();
        }
    }

    /* Now add all of the edges */
    for (const auto &edge : edges) {
        addEdge(edge.u, edge.v, edge.weight);
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateUndirected(NodeT u, size_t size) {
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

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateDirected(NodeT u, size_t outSize, size_t inSize) {
    preallocateDirectedOutEdges(u, outSize);
    preallocateDirectedInEdges(u, inSize);
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateDirectedOutEdges(NodeT u, size_t outSize) {
    assert(directed);
    assert(exists[u]);
    outEdges[u].reserve(outSize);

    if (weighted) {
        outEdgeWeights[u].reserve(outSize);
    }
    if (edgesIndexed) {
        outEdges[u].reserve(outSize);
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateDirectedInEdges(NodeT u, size_t inSize) {
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
/** PRIVATE HELPERS **/

template <class NodeT, class EdgeWeightT>
index AdjListGraph<NodeT, EdgeWeightT>::indexInInEdgeArray(NodeT v, NodeT u) const {
    if (!directed) {
        return indexInOutEdgeArray(v, u);
    }
    for (index i = 0; i < inEdges[v].size(); ++i) {
        NodeT x = inEdges[v][i];
        if (x == u) {
            return i;
        }
    }
    return none;
}

template <class NodeT, class EdgeWeightT>
index AdjListGraph<NodeT, EdgeWeightT>::indexInOutEdgeArray(NodeT u, NodeT v) const {
    for (index i = 0; i < outEdges[u].size(); ++i) {
        NodeT x = outEdges[u][i];
        if (x == v) {
            return i;
        }
    }
    return none;
}

/** EDGE IDS **/

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::indexEdges(bool force) {
    if (edgesIndexed && !force)
        return;

    omega = 0; // reset edge ids (for re-indexing)

    outEdgeIds.clear(); // reset ids vector (for re-indexing)
    outEdgeIds.resize(outEdges.size());
    forNodes([&](NodeT u) { outEdgeIds[u].resize(outEdges[u].size(), nullEdgeId); });

    if (directed) {
        inEdgeIds.resize(inEdges.size());
        forNodes([&](NodeT u) { inEdgeIds[u].resize(inEdges[u].size(), nullEdgeId); });
    }

    // assign edge ids for edges in one direction
    forNodes([&](NodeT u) {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            NodeT v = outEdges[u][i];
            if (v != nullNodeId && (directed || (u >= v))) {
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
        balancedParallelForNodes([&](NodeT u) {
            for (index i = 0; i < outEdges[u].size(); ++i) {
                NodeT v = outEdges[u][i];
                if (v != nullNodeId && outEdgeIds[u][i] == nullEdgeId) {
                    index j = indexInOutEdgeArray(v, u);
                    outEdgeIds[u][i] = outEdgeIds[v][j];
                }
            }
        });
    } else {
        balancedParallelForNodes([&](NodeT u) {
            for (index i = 0; i < inEdges[u].size(); ++i) {
                NodeT v = inEdges[u][i];
                if (v != nullNodeId) {
                    index j = indexInOutEdgeArray(v, u);
                    inEdgeIds[u][i] = outEdgeIds[v][j];
                }
            }
        });
    }

    edgesIndexed = true; // remember that edges have been indexed so that addEdge
                         // needs to create edge ids
}

template <class NodeT, class EdgeWeightT>
edgeid AdjListGraph<NodeT, EdgeWeightT>::edgeId(NodeT u, NodeT v) const {
    if (!edgesIndexed) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    const index i = indexInOutEdgeArray(u, v);

    if (i == none) {
        throw std::runtime_error("Edge does not exist");
    }
    return outEdgeIds[u][i];
}

/** GRAPH INFORMATION **/

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::shrinkToFit() {
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

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::compactEdges() {
    this->parallelForNodes([&](NodeT u) {
        if (degreeOut(u) == 0) {
            outEdges[u].clear();
            if (weighted)
                outEdgeWeights[u].clear();
            if (edgesIndexed)
                outEdgeIds[u].clear();
        } else {
            for (index i = 0; i < outEdges[u].size(); ++i) {
                while (i < outEdges[u].size() && outEdges[u][i] == nullNodeId) {
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
                    while (i < inEdges[u].size() && inEdges[u][i] == nullNodeId) {
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

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::sortEdges() {
    std::vector<std::vector<NodeT>> targetAdjacencies(upperNodeIdBound());
    std::vector<std::vector<EdgeWeightT>> targetWeight;
    std::vector<std::vector<edgeid>> targetEdgeIds;

    if (isWeighted()) {
        targetWeight.resize(upperNodeIdBound());
        forNodes([&](NodeT u) { targetWeight[u].reserve(degree(u)); });
    }
    if (hasEdgeIds()) {
        targetEdgeIds.resize(upperNodeIdBound());
        forNodes([&](NodeT u) { targetEdgeIds[u].reserve(degree(u)); });
    }

    forNodes([&](NodeT u) { targetAdjacencies[u].reserve(degree(u)); });

    auto assignToTarget = [&](NodeT u, NodeT v, EdgeWeightT w, edgeid eid) {
        targetAdjacencies[v].push_back(u);
        if (isWeighted()) {
            targetWeight[v].push_back(w);
        }
        if (hasEdgeIds()) {
            targetEdgeIds[v].push_back(eid);
        }
    };

    forNodes([&](NodeT u) { forInEdgesOf(u, assignToTarget); });

    outEdges.swap(targetAdjacencies);
    outEdgeWeights.swap(targetWeight);
    outEdgeIds.swap(targetEdgeIds);

    if (isDirected()) {
        inEdges.swap(targetAdjacencies);
        inEdgeWeights.swap(targetWeight);
        inEdgeIds.swap(targetEdgeIds);

        forNodes([&](NodeT u) {
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

        forNodes([&](NodeT u) { forEdgesOf(u, assignToTarget); });

        inEdges.swap(targetAdjacencies);
        inEdgeWeights.swap(targetWeight);
        inEdgeIds.swap(targetEdgeIds);
    }
}

template <class NodeT, class EdgeWeightT>
EdgeWeightT
AdjListGraph<NodeT, EdgeWeightT>::computeWeightedDegree(NodeT u, bool inDegree,
                                                        bool countSelfLoopsTwice) const {
    if (weighted) {
        EdgeWeightT sum{0};
        auto sumWeights = [&](NodeT v, EdgeWeightT w) {
            sum += (countSelfLoopsTwice && u == v) ? EdgeWeightT{2} * w : w;
        };
        if (inDegree) {
            forInNeighborsOf(u, sumWeights);
        } else {
            forNeighborsOf(u, sumWeights);
        }
        return sum;
    }

    count sum = inDegree ? degreeIn(u) : degreeOut(u);
    auto countSelfLoops = [&](NodeT v) { sum += (u == v); };

    if (countSelfLoopsTwice && numberOfSelfLoops()) {
        if (inDegree) {
            forInNeighborsOf(u, countSelfLoops);
        } else {
            forNeighborsOf(u, countSelfLoops);
        }
    }

    return static_cast<EdgeWeightT>(sum);
}

/** NODE MODIFIERS **/

template <class NodeT, class EdgeWeightT>
NodeT AdjListGraph<NodeT, EdgeWeightT>::addNode() {
    NodeT v = z; // node gets maximum id
    ++z;         // increment node range
    ++n;         // increment node count

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

template <class NodeT, class EdgeWeightT>
NodeT AdjListGraph<NodeT, EdgeWeightT>::addNodes(count numberOfNewNodes) {
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

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeNode(NodeT v) {
    assert(v < z);
    assert(exists[v]);

    // Remove all outgoing and ingoing edges
    while (!outEdges[v].empty())
        removeEdge(v, outEdges[v].front());
    if (isDirected())
        while (!inEdges[v].empty())
            removeEdge(inEdges[v].front(), v);

    // Make the attributes of this node invalid
    auto &theMap = nodeAttributeMap.attrMap;
    for (auto it = theMap.begin(); it != theMap.end(); ++it) {
        auto attributeStorageBase = it->second.get();
        attributeStorageBase->invalidate(v);
    }

    exists[v] = false;
    --n;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::restoreNode(NodeT v) {
    assert(v < z);
    assert(!exists[v]);

    exists[v] = true;
    ++n;
}

/** NODE PROPERTIES **/

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::weightedDegree(NodeT u,
                                                             bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, false, countSelfLoopsTwice);
}

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::weightedDegreeIn(NodeT u,
                                                               bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, true, countSelfLoopsTwice);
}

/** EDGE MODIFIERS **/

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addEdge(NodeT u, NodeT v, EdgeWeightT ew,
                                               bool checkForMultiEdges) {
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

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addPartialEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew,
                                                      uint64_t index, bool checkForMultiEdges) {
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

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addPartialOutEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew,
                                                         uint64_t index, bool checkForMultiEdges) {
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

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addPartialInEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew,
                                                        uint64_t index, bool checkForMultiEdges) {
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

template <class NodeT, typename T>
void erase(NodeT u, index idx, std::vector<std::vector<T>> &vec) {
    vec[u][idx] = vec[u].back();
    vec[u].pop_back();
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeEdge(NodeT u, NodeT v) {
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
    --m; // decrease number of edges
    if (isLoop)
        --storedNumberOfSelfLoops;

    // remove edge for source node
    erase(u, vi, outEdges);
    if (weighted) {
        erase(u, vi, outEdgeWeights);
    }
    if (edgesIndexed) {
        erase(u, vi, outEdgeIds);
        // Make the attributes of this edge invalid
        auto &theMap = edgeAttributeMap.attrMap;
        for (auto it = theMap.begin(); it != theMap.end(); ++it) {
            auto attributeStorageBase = it->second.get();
            attributeStorageBase->invalidate(deletedID);
        }
    }
    if (!directed && !isLoop) {
        // also remove edge for target node
        erase(v, ui, outEdges);
        if (weighted) {
            erase(v, ui, outEdgeWeights);
        }
        if (edgesIndexed) {
            erase(v, ui, outEdgeIds);
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
                auto &theMap = edgeAttributeMap.attrMap;
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
        balancedParallelForNodes([&](NodeT w) {
            for (index i = 0; i < outEdges[w].size(); ++i) {
                auto curID = outEdgeIds[w][i];
                if (curID > deletedID) {
                    --outEdgeIds[w][i];
                }
            }
        });
        // use erase to remove data entry at index `deletedID` and compact the data vector again
        auto &theMap = edgeAttributeMap.attrMap;
        for (auto it = theMap.begin(); it != theMap.end(); ++it) {
            auto attributeStorageBase = it->second.get();
            attributeStorageBase->erase(deletedID);
        }
    }
    if (directed) {
        assert(ui != none);

        erase(v, ui, inEdges);
        if (weighted) {
            erase(v, ui, inEdgeWeights);
        }
        if (edgesIndexed) {
            erase(v, ui, inEdgeIds);
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
            balancedParallelForNodes([&](NodeT w) {
                for (index i = 0; i < inEdges[w].size(); ++i) {
                    NodeT vv = inEdges[w][i];
                    if (vv != none) {
                        index j = indexInOutEdgeArray(vv, w);
                        inEdgeIds[w][i] = outEdgeIds[vv][j];
                    }
                }
            });
        }
    }
    if (maintainCompactEdges) {
        --omega; // decrease upperBound of edges
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeAllEdges() {
    parallelForNodes([&](const NodeT u) {
        removePartialOutEdges(unsafe, u);
        if (isDirected()) {
            removePartialInEdges(unsafe, u);
        }
    });

    m = 0;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeSelfLoops() {
    parallelForNodes([&](const NodeT u) {
        auto isSelfLoop = [u](const NodeT v) { return u == v; };
        removeAdjacentEdges(u, isSelfLoop);
        if (isDirected()) {
            removeAdjacentEdges(u, isSelfLoop, true);
        }
    });

    m -= storedNumberOfSelfLoops;
    storedNumberOfSelfLoops = 0;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeMultiEdges() {
    count removedEdges = 0;
    count removedSelfLoops = 0;
    std::unordered_set<NodeT> nodes;

    forNodes([&](const NodeT u) {
        nodes.reserve(degree(u));
        auto isMultiedge = [&nodes](const NodeT v) { return !nodes.insert(v).second; };
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

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::swapEdge(NodeT s1, NodeT t1, NodeT s2, NodeT t2) {
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

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::hasEdge(NodeT u, NodeT v) const noexcept {
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
}

/** EDGE ATTRIBUTES **/

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::weight(NodeT u, NodeT v) const {
    index vi = indexInOutEdgeArray(u, v);
    if (vi == none) {
        return nullWeight;
    } else {
        return weighted ? outEdgeWeights[u][vi] : defaultEdgeWeightT;
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::setWeight(NodeT u, NodeT v, EdgeWeightT ew) {
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

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::increaseWeight(NodeT u, NodeT v, EdgeWeightT ew) {
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

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::setWeightAtIthNeighbor(Unsafe, NodeT u, index i,
                                                              EdgeWeightT ew) {
    outEdgeWeights[u][i] = ew;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::setWeightAtIthInNeighbor(Unsafe, NodeT u, index i,
                                                                EdgeWeightT ew) {
    inEdgeWeights[u][i] = ew;
}

/** SUMS **/

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::totalEdgeWeight() const noexcept {
    if (weighted)
        return parallelSumForEdges([](NodeT, NodeT, EdgeWeightT ew) { return ew; });
    return numberOfEdges() * defaultEdgeWeightT;
}

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::checkConsistency() const {
    // check for multi-edges
    std::vector<NodeT> lastSeen(z, nullNodeId);
    bool noMultiEdges = true;
    auto noMultiEdgesDetected = [&noMultiEdges]() { return noMultiEdges; };
    forNodesWhile(noMultiEdgesDetected, [&](NodeT v) {
        forNeighborsOf(v, [&](NodeT u) {
            if (lastSeen[u] == v) {
                noMultiEdges = false;
                DEBUG("Multiedge found between ", u, " and ", v, "!");
            }
            lastSeen[u] = v;
        });
    });

    bool correctNodeUpperbound = (z == outEdges.size()) && ((directed ? z : 0) == inEdges.size())
                                 && ((weighted ? z : 0) == outEdgeWeights.size())
                                 && ((weighted && directed ? z : 0) == inEdgeWeights.size())
                                 && ((edgesIndexed ? z : 0) == outEdgeIds.size())
                                 && ((edgesIndexed && directed ? z : 0) == inEdgeIds.size());

    if (!correctNodeUpperbound)
        DEBUG("Saved node upper bound doesn't actually match the actual node upper bound!");

    count NumberOfOutEdges = 0;
    count NumberOfOutEdgeWeights = 0;
    count NumberOfOutEdgeIds = 0;
    for (index i = 0; i < outEdges.size(); ++i) {
        NumberOfOutEdges += outEdges[i].size();
    }
    if (weighted)
        for (index i = 0; i < outEdgeWeights.size(); ++i) {
            NumberOfOutEdgeWeights += outEdgeWeights[i].size();
        }
    if (edgesIndexed)
        for (index i = 0; i < outEdgeIds.size(); ++i) {
            NumberOfOutEdgeIds += outEdgeIds[i].size();
        }

    count NumberOfInEdges = 0;
    count NumberOfInEdgeWeights = 0;
    count NumberOfInEdgeIds = 0;
    if (directed) {
        for (index i = 0; i < inEdges.size(); ++i) {
            NumberOfInEdges += inEdges[i].size();
        }
        if (weighted)
            for (index i = 0; i < inEdgeWeights.size(); ++i) {
                NumberOfInEdgeWeights += inEdgeWeights[i].size();
            }
        if (edgesIndexed)
            for (index i = 0; i < inEdgeIds.size(); ++i) {
                NumberOfInEdgeIds += inEdgeIds[i].size();
            }
    }

    if (!directed) {
        NumberOfOutEdges = (NumberOfOutEdges + storedNumberOfSelfLoops) / 2;
        if (weighted)
            NumberOfOutEdgeWeights = (NumberOfOutEdgeWeights + storedNumberOfSelfLoops) / 2;
        if (edgesIndexed)
            NumberOfOutEdgeIds = (NumberOfOutEdgeIds + storedNumberOfSelfLoops) / 2;
    }

    bool correctNumberOfEdges = (m == NumberOfOutEdges) && ((directed ? m : 0) == NumberOfInEdges)
                                && ((weighted ? m : 0) == NumberOfOutEdgeWeights)
                                && ((weighted && directed ? m : 0) == NumberOfInEdgeWeights)
                                && ((edgesIndexed ? m : 0) == NumberOfOutEdgeIds)
                                && ((edgesIndexed && directed ? m : 0) == NumberOfInEdgeIds);

    if (!correctNumberOfEdges)
        DEBUG("Saved number of edges is incorrect!");

    return noMultiEdges && correctNodeUpperbound && correctNumberOfEdges;
}

template <class NodeT, class EdgeWeightT>
AdjListGraph<NodeT, EdgeWeightT> AdjListGraph<NodeT, EdgeWeightT>::fromCSR(
    std::span<const index> rowIdxView, std::span<const index> columnIdxView,
    std::span<const double> nonZerosView, bool directed, bool isWeighted) {
    // In CSR format rowIdxView has nRows + 1 entries, so the number of rows is
    // derived directly from its length.
    if (rowIdxView.empty()) {
        throw std::invalid_argument("rowIdxView must have at least one entry");
    }
    count nRows = static_cast<count>(rowIdxView.size() - 1);

    // In CSR format the number of stored entries (nnz) is the last row pointer.
    size_t nnz = static_cast<size_t>(rowIdxView[nRows]);

    if (columnIdxView.size() != nnz) {
        throw std::invalid_argument("columnIdxView size does not match rowIdxView[nRows]");
    }
    if (isWeighted && nonZerosView.size() != nnz) {
        throw std::invalid_argument("nonZerosView size does not match rowIdxView[nRows]");
    }

    // Create the Graph with reserved structures
    AdjListGraph graph(nRows, isWeighted, directed);

    // Number of edges equals number of stored column indices
    count edgeCount = static_cast<count>(nnz);

    // Count self-loops so that numberOfSelfLoops() stays consistent with the
    // graph that addEdge() would have produced from the same entries.
    count selfLoops = 0;
    for (index i = 0; i < nRows; ++i) {
        assert(rowIdxView[i + 1] >= rowIdxView[i] && "rowIdxView must be non-decreasing");
        for (index k = rowIdxView[i]; k < rowIdxView[i + 1]; ++k) {
            if (columnIdxView[k] == i)
                ++selfLoops;
        }
    }
    graph.storedNumberOfSelfLoops = selfLoops;

    graph.outEdges.resize(nRows);
    if (isWeighted)
        graph.outEdgeWeights.resize(nRows);

    if (directed) {
        // Directed: row i lists the out-neighbors of node i directly.
        for (index i = 0; i < nRows; ++i) {
            size_t rowSize = static_cast<size_t>(rowIdxView[i + 1] - rowIdxView[i]);
            graph.outEdges[i].resize(rowSize);
            if (isWeighted)
                graph.outEdgeWeights[i].resize(rowSize);
        }

#pragma omp parallel for schedule(guided)
        for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
            index rowStart = rowIdxView[i];
            index rowEnd = rowIdxView[i + 1];
            if (rowEnd > rowStart) {
                std::copy(columnIdxView.data() + rowStart, columnIdxView.data() + rowEnd,
                          graph.outEdges[i].begin());
                if (isWeighted) {
                    std::copy(nonZerosView.data() + rowStart, nonZerosView.data() + rowEnd,
                              graph.outEdgeWeights[i].begin());
                }
            }
        }

        // Build inEdges (and inEdgeWeights). We compute in-degrees first,
        // preallocate, then fill the arrays. Filling uses atomic position
        // counters so it can be parallelized safely.
        graph.inEdges.resize(nRows);
        if (isWeighted) {
            graph.inEdgeWeights.resize(nRows);
        }

        std::vector<size_t> inDeg(nRows, 0);
        for (size_t k = 0; k < nnz; ++k) {
            index j = columnIdxView[k];
            ++inDeg[j];
        }

        for (index j = 0; j < nRows; ++j) {
            graph.inEdges[j].resize(inDeg[j]);
            if (isWeighted)
                graph.inEdgeWeights[j].resize(inDeg[j]);
        }

        // position counters for each target row; initialize to zero
        std::vector<std::atomic<size_t>> inPos(nRows);
        for (index j = 0; j < nRows; ++j)
            inPos[j].store(0);

#pragma omp parallel for schedule(guided)
        for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
            index rowStart = rowIdxView[i];
            index rowEnd = rowIdxView[i + 1];
            for (index k = rowStart; k < rowEnd; ++k) {
                index j = columnIdxView[k];
                size_t pos = inPos[j].fetch_add(1);
                graph.inEdges[j][pos] = static_cast<NodeT>(i);
                if (isWeighted)
                    graph.inEdgeWeights[j][pos] = nonZerosView[k];
            }
        }
    } else {
        // Undirected: every stored nonzero (i, j) is an undirected edge and is
        // stored from both endpoints (self-loops are stored only once), exactly
        // as addEdge() does.
        //
        // First determine the final degree of each node: every entry adds a
        // slot to its row owner, and every non-self-loop entry adds a mirrored
        // slot to its neighbor.
        std::vector<size_t> deg(nRows, 0);
        for (index i = 0; i < nRows; ++i) {
            for (index k = rowIdxView[i]; k < rowIdxView[i + 1]; ++k) {
                index j = columnIdxView[k];
                ++deg[i];
                if (j != i)
                    ++deg[j];
            }
        }

        for (index i = 0; i < nRows; ++i) {
            graph.outEdges[i].resize(deg[i]);
            if (isWeighted)
                graph.outEdgeWeights[i].resize(deg[i]);
        }

        // Per-node position counters. A node's slots are written both by the
        // thread owning its row and by threads mirroring edges into it, so the
        // counters must be atomic.
        std::vector<std::atomic<size_t>> pos(nRows);
        for (index i = 0; i < nRows; ++i)
            pos[i].store(0);

#pragma omp parallel for schedule(guided)
        for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
            for (index k = rowIdxView[i]; k < rowIdxView[i + 1]; ++k) {
                index j = columnIdxView[k];
                bool mirror = j != static_cast<index>(i);

                size_t p = pos[i].fetch_add(1);
                graph.outEdges[i][p] = static_cast<NodeT>(j);

                size_t q = 0;
                if (mirror) {
                    q = pos[j].fetch_add(1);
                    graph.outEdges[j][q] = static_cast<NodeT>(i);
                }

                if (isWeighted) {
                    EdgeWeightT w = static_cast<EdgeWeightT>(nonZerosView[k]);
                    graph.outEdgeWeights[i][p] = w;
                    if (mirror)
                        graph.outEdgeWeights[j][q] = w;
                }
            }
        }
    }

    graph.m = edgeCount;

    return graph;
}

template <class NodeT, class EdgeWeightT>
AdjListGraph<NodeT, EdgeWeightT>
AdjListGraph<NodeT, EdgeWeightT>::_fromCSRRaw(const index *rowIdxPtr, std::size_t rowIdxSize,
                                              const index *columnIdxPtr, std::size_t columnIdxSize,
                                              const double *nonZerosPtr, std::size_t nonZerosSize,
                                              bool directed, bool isWeighted) {
    // This helper only adapts raw pointers to std::span for callers (e.g. Cython)
    // that cannot construct a span directly; the real work happens in fromCSR().
    if (rowIdxPtr == nullptr || columnIdxPtr == nullptr) {
        throw std::invalid_argument("null CSR pointer passed to _fromCSRRaw");
    }
    if (isWeighted && nonZerosPtr == nullptr) {
        throw std::invalid_argument(
            "null nonZeros pointer passed to _fromCSRRaw for weighted graph");
    }
    return fromCSR(std::span<const index>(rowIdxPtr, rowIdxSize),
                   std::span<const index>(columnIdxPtr, columnIdxSize),
                   std::span<const double>(nonZerosPtr, nonZerosSize), directed, isWeighted);
}

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_ADJ_LIST_GRAPH_IMPL_HPP_
