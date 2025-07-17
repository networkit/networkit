#ifndef NETWORKIT_GRAPH_GRAPHW_HPP_
#define NETWORKIT_GRAPH_GRAPHW_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * A writable graph that extends Graph with mutation operations.
 * This class provides all read operations from Graph plus write operations
 * like addNode, addEdge, removeNode, removeEdge, etc.
 */
class GraphW final : public Graph {
public:
    /**
     * Create a graph of @a n nodes. The graph has assignable edge weights if @a
     * weighted is set to <code>true</code>. If @a weighted is set to
     * <code>false</code> each edge has edge weight 1.0 and any other weight
     * assignment will be ignored.
     * @param n Number of nodes.
     * @param weighted If set to <code>true</code>, the graph has edge weights.
     * @param directed If set to @c true, the graph will be directed.
     * @param edgesIndexed If set to @c true, the graph will have indexed edges.
     */
    GraphW(count n = 0, bool weighted = false, bool directed = false, bool edgesIndexed = false)
        : Graph(n, weighted, directed, edgesIndexed) {}

    /**
     * Generate a weighted graph from a list of edges. (Useful for small
     * graphs in unit tests that you do not want to read from a file.)
     *
     * @param[in] edges list of weighted edges
     */
    GraphW(std::initializer_list<WeightedEdge> edges);

    /**
     * Create a graph as copy of @a other.
     * @param other The graph to copy.
     */
    GraphW(const GraphW &other) : Graph(other) {}

    /**
     * Create a graph as copy of @a other.
     * @param other The graph to copy.
     */
    GraphW(const Graph &other) : Graph(other) {}

    /**
     * Create a graph as copy of @a other with modified properties.
     * @param other The graph to copy.
     * @param weighted If set to true, the graph has edge weights.
     * @param directed If set to true, the graph will be directed.
     * @param edgesIndexed If set to true, the graph will have indexed edges.
     */
    template <class EdgeMerger = std::plus<edgeweight>>
    GraphW(const Graph &other, bool weighted, bool directed, bool edgesIndexed = false,
           EdgeMerger edgeMerger = std::plus<edgeweight>())
        : Graph(other, weighted, directed, edgesIndexed, edgeMerger) {}

    /** move constructor */
    GraphW(GraphW &&other) noexcept : Graph(std::move(other)) {}

    /** move constructor */
    GraphW(Graph &&other) noexcept : Graph(std::move(other)) {}

    /** Default destructor */
    ~GraphW() = default;

    /** move assignment operator */
    GraphW &operator=(GraphW &&other) noexcept {
        Graph::operator=(std::move(other));
        return *this;
    }

    /** move assignment operator */
    GraphW &operator=(Graph &&other) noexcept {
        Graph::operator=(std::move(other));
        return *this;
    }

    /** copy assignment operator */
    GraphW &operator=(const GraphW &other) {
        Graph::operator=(other);
        return *this;
    }

    /** copy assignment operator */
    GraphW &operator=(const Graph &other) {
        Graph::operator=(other);
        return *this;
    }

    /** EDGE IDS **/

    /**
     * Initially assign integer edge identifiers.
     *
     * @param force Force re-indexing of edges even if they have already been
     * indexed
     */
    void indexEdges(bool force = false);

    /** GRAPH INFORMATION **/

    /**
     * Try to save some memory by shrinking internal data structures of the
     * graph. Only run this once you finished editing the graph. Otherwise it
     * will cause unnecessary reallocation of memory.
     */
    void shrinkToFit();

    /**
     * DEPRECATED: this function will no longer be supported in later releases.
     * Compacts the adjacency arrays by re-using no longer needed slots from
     * deleted edges.
     */
    void TLX_DEPRECATED(compactEdges());

    /**
     * Sorts the outgoing neighbors of a given node according to a user-defined comparison function.
     *
     * @param u The node whose outgoing neighbors will be sorted.
     * @param lambda A binary predicate used to compare two neighbors. The predicate should
     *               take two nodes as arguments and return true if the first node should
     *               precede the second in the sorted order.
     */
    template <typename Lambda>
    void sortNeighbors(node u, Lambda lambda) {
        assert(hasNode(u));
        std::sort(outEdges[u].begin(), outEdges[u].end(), lambda);
        if (isDirected()) {
            std::sort(inEdges[u].begin(), inEdges[u].end(), lambda);
        }
    }

    /**
     * Sorts the adjacency arrays by node id. While the running time is linear
     * this temporarily duplicates the memory.
     */
    void sortEdges();

    /**
     * Sorts the adjacency arrays by a custom criterion.
     *
     * @param lambda Lambda function used to sort the edges. It takes two WeightedEdge
     * e1 and e2 as input parameters, returns true if e1 < e2, false otherwise.
     */
    template <class Lambda>
    void sortEdges(Lambda lambda) {
        parallelForNodes([&](node u) {
            if (isWeighted()) {
                std::vector<WeightedEdge> edges;
                forNeighborsOf(u, [&](node v, edgeweight w) {
                    edges.emplace_back(u, v, w);
                });
                std::sort(edges.begin(), edges.end(), lambda);

                removePartialOutEdges(unsafe, u);
                for (const auto& edge : edges) {
                    addPartialOutEdge(unsafe, edge.u, edge.v, edge.weight);
                }
            } else {
                std::vector<node> neighbors(outEdges[u]);
                std::sort(neighbors.begin(), neighbors.end(), [&](node v1, node v2) {
                    return lambda(WeightedEdge(u, v1, defaultEdgeWeight),
                                 WeightedEdge(u, v2, defaultEdgeWeight));
                });

                removePartialOutEdges(unsafe, u);
                for (node v : neighbors) {
                    addPartialOutEdge(unsafe, u, v);
                }
            }
        });
    }

    /**
     * Set edge count of the graph to edges.
     * @param edges the edge count of a graph
     */
    void setEdgeCount(Unsafe, count edges) { m = edges; }

    /**
     * Set upper bound of edge count.
     *
     * @param newBound New upper edge id bound.
     */
    void setUpperEdgeIdBound(Unsafe, edgeid newBound) { omega = newBound; }

    /**
     * Set the number of self-loops.
     *
     * @param loops New number of self-loops.
     */
    void setNumberOfSelfLoops(Unsafe, count loops) { storedNumberOfSelfLoops = loops; }

    /* NODE MODIFIERS */

    /**
     * Add a new node to the graph and return it.
     * @return The new node.
     */
    node addNode();

    /**
     * Add numberOfNewNodes new nodes.
     * @param  numberOfNewNodes Number of new nodes.
     * @return The index of the last node added.
     */
    node addNodes(count numberOfNewNodes);

    /**
     * Remove a node @a v and all incident edges from the graph.
     *
     * Incoming as well as outgoing edges will be removed.
     *
     * @param v Node.
     */
    void removeNode(node v);

    /**
     * Removes out-going edges from node @u. If the graph is weighted and/or has edge ids, weights
     * and/or edge ids will also be removed.
     *
     * @param u Node.
     */
    void removePartialOutEdges(Unsafe, node u) {
        assert(hasNode(u));
        outEdges[u].clear();
        if (isWeighted()) {
            outEdgeWeights[u].clear();
        }
        if (hasEdgeIds()) {
            outEdgeIds[u].clear();
        }
    }

    /**
     * Removes in-going edges to node @u. If the graph is weighted and/or has edge ids, weights
     * and/or edge ids will also be removed.
     *
     * @param u Node.
     */
    void removePartialInEdges(Unsafe, node u) {
        assert(hasNode(u));
        inEdges[u].clear();
        if (isWeighted()) {
            inEdgeWeights[u].clear();
        }
        if (hasEdgeIds()) {
            inEdgeIds[u].clear();
        }
    }

    /**
     * Restores a previously deleted node @a v with its previous id in the
     * graph.
     *
     * @param v Node.
     *
     */
    void restoreNode(node v);

    /* EDGE MODIFIERS */

    /**
     * Insert an edge between the nodes @a u and @a v. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param checkMultiEdge If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addEdge(node u, node v, edgeweight ew = defaultEdgeWeight, bool checkMultiEdge = false);

    /**
     * Insert an edge between the nodes @a u and @a v. Unlike the addEdge function, this function
     * does not add any information to v. If the graph is weighted you can optionally set a
     * weight for this edge. The default weight is 1.0. Note: Multi-edges are not supported and will
     * NOT be handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param index Optional edge index.
     * @param checkForMultiEdges If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight,
                        uint64_t index = 0, bool checkForMultiEdges = false);

    /**
     * Insert an in edge between the nodes @a u and @a v in a directed graph. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param index Optional edge index.
     * @param checkForMultiEdges If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialInEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight,
                          uint64_t index = 0, bool checkForMultiEdges = false);

    /**
     * Insert an out edge between the nodes @a u and @a v in a directed graph. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param index Optional edge index.
     * @param checkForMultiEdges If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialOutEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight,
                           uint64_t index = 0, bool checkForMultiEdges = false);

    /**
     * If set to true, the ingoing and outgoing adjacency vectors will
     * automatically be updated to maintain a sorting (if it existed before) by performing up to n-1
     * swaps. If the user plans to remove multiple edges, better set it to false and call
     * sortEdges() afterwards to avoid redundant swaps. Default = true.
     */
    void setKeepEdgesSorted(bool sorted = true) { maintainSortedEdges = sorted; }

    /**
     * If set to true, the EdgeIDs will automatically be adjusted,
     * so that no gaps in between IDs exist. If the user plans to remove multiple edges, better set
     * it to false and call indexEdges(force=true) afterwards to avoid redundant re-indexing.
     * Default = true.
     */
    void setMaintainCompactEdges(bool compact = true) { maintainCompactEdges = compact; }

    /**
     * Removes the undirected edge {@a u,@a v}.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     */
    void removeEdge(node u, node v);

    /**
     * Removes all the edges in the graph.
     */
    void removeAllEdges();

    /**
     * Removes edges adjacent to a node according to a specific criterion.
     *
     * @param u The node whose adjacent edges shall be removed.
     * @param condition A function that takes a node as an input and returns a
     * bool. If true the edge (u, v) is removed.
     * @param edgesIn Whether in-going or out-going edges shall be removed.
     * @return std::pair<count, count> The number of removed edges (first) and the number of removed
     * self-loops (second).
     */
    template <typename Condition>
    std::pair<count, count> removeAdjacentEdges(node u, Condition condition, bool edgesIn = false) {
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

    /**
     * Removes all self-loops in the graph.
     */
    void removeSelfLoops();

    /**
     * Removes all multi-edges in the graph.
     */
    void removeMultiEdges();

    /**
     * Changes the edges {@a s1, @a t1} into {@a s1, @a t2} and the edge {@a
     * s2,
     * @a t2} into {@a s2, @a t1}.
     *
     * If there are edge weights or edge ids, they are preserved. Note that no
     * check is performed if the swap is actually possible, i.e. does not
     * generate duplicate edges.
     *
     * @param s1 The first source
     * @param t1 The first target
     * @param s2 The second source
     * @param t2 The second target
     */
    void swapEdge(node s1, node t1, node s2, node t2);

    /**
     * Set edge weight of edge {@a u,@a v}. BEWARE: Running time is \Theta(deg(u))!
     *
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew New edge weight.
     */
    void setWeight(node u, node v, edgeweight ew);

    /**
     * Set edge weight of the @a i-th outgoing edge of node @a u. BEWARE: Running time is constant.
     *
     * @param u Endpoint of edge.
     * @param i Index of the outgoing edge.
     * @param ew New edge weight.
     */
    void setWeightAtIthNeighbor(Unsafe, node u, index i, edgeweight ew);

    /**
     * Set edge weight of the @a i-th incoming edge of node @a u. BEWARE: Running time is constant.
     *
     * @param u Endpoint of edge.
     * @param i Index of the incoming edge.
     * @param ew New edge weight.
     */
    void setWeightAtIthInNeighbor(Unsafe, node u, index i, edgeweight ew);

    /**
     * Increase edge weight of edge {@a u,@a v} by @a ew. BEWARE: Running time is \Theta(deg(u))!
     *
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Edge weight increase.
     */
    void increaseWeight(node u, node v, edgeweight ew);

    /**
     * Reserves memory in the node's edge containers for undirected graphs.
     *
     * @param u the node memory should be reserved for
     * @param size the amount of memory to reserve
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateUndirected(node u, size_t size);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param inSize the amount of memory to reserve for in edges
     * @param outSize the amount of memory to reserve for out edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirected(node u, size_t outSize, size_t inSize);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param outSize the amount of memory to reserve for out edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirectedOutEdges(node u, size_t outSize);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param inSize the amount of memory to reserve for in edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirectedInEdges(node u, size_t inSize);
};

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_GRAPHW_HPP_
