#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H

#include <unordered_map>
#include <networkit/graph/Graph.hpp>
#include <tlx/unused.hpp>

namespace NetworKit {
namespace GraphTools {

/**
 * Computes a graph with the same structure but with continuous node ids.
 * @param  graph     The graph to be compacted.
 * @param  nodeIdMap The map providing the information about the node ids.
 * @return           Returns a compacted Graph.
 */
Graph getCompactedGraph(const Graph& graph, const std::unordered_map<node,node>& nodeIdMap);

/**
 * Computes a map of node ids.
 * @param	graph	The graph of which the node id map is wanted.
 * @return			Returns the node id map.
 */
std::unordered_map<node,node> getContinuousNodeIds(const Graph& graph);

/**
 * Computes a map of random node ids
 * @param	graph	The graph of which the node id map is wanted.
 * @return		Returns the node id map.
 */
std::unordered_map<node, node> getRandomContinuousNodeIds(const Graph& graph);


/**
 * Inverts a given mapping of node ids from a graph with deleted nodes to continuous node ids.
 * @param 	nodeIdMap	The mapping from node ids with gaps to continuous node ids (i.e. from @getContinuousNodeIds)
 * @param 	G 			The compacted graph (currently only needed for the upper node id bound)
 * @return 				A vector of nodes id where the index is the node id of the compacted graph and the value is the node id of the noncontinuous graph.
 */
std::vector<node> invertContinuousNodeIds(const std::unordered_map<node,node>& nodeIdMap, const Graph& G);

/**
 * Constructs a new graph that has the same node ids as before it was compacted.
 * @param  invertedIdMap The node id mapping from continuous node ids to noncontinuous node ids.
 * @param  G             The compacted graph.
 * @return               The original graph.
 */
Graph restoreGraph(const std::vector<node>& invertedIdMap, const Graph& G);


/**
 * Rename nodes in a graph using a callback which translates each old id to a new one.
 * For each node u in input graph, oldIdToNew(u) < numNodes.
 *
 * @param graph Input graph.
 * @param numNodes    Number of nodes in the output graph.
 * @param oldIdToNew  Translate old id to new ones. Must be thread-safe
 * @param skipNode    Skip all nodes (and incident edges) for old node
 *                    ids u where deleteNode(u) == true, Must be thread-safe
 * @param preallocate Preallocates memory before adding neighbors
 *                    (Preallocation does not account for deleted nodes
 *                    and hence may need more memory)
 *
 * @node preallocate is currently not implemented
 */
template <typename UnaryIdMapper, typename SkipEdgePredicate>
Graph getRemappedGraph(const Graph& graph, count numNodes,
    UnaryIdMapper&& oldIdToNew, SkipEdgePredicate&& skipNode, bool preallocate = true)
{
    tlx::unused(preallocate); // TODO: Add perallocate as soon as Graph supports it

#ifndef NDEBUG
    graph.forNodes([&] (node u) {
        assert(skipNode(u) || oldIdToNew(u) < numNodes);
    });
#endif

    const auto directed = graph.isDirected();
    Graph Gnew(numNodes, graph.isWeighted(), directed);

    graph.forNodes([&](const node u) { // TODO: Make parallel when graph support addHalfEdge
        if (skipNode(u))
            return;

        const node mapped_u = oldIdToNew(u);
        graph.forNeighborsOf(u, [&](node, node v, edgeweight ew) {
            if (!directed && v < u)
                return;
            if (skipNode(v))
                return;

            const node mapped_v = oldIdToNew(v);
            Gnew.addEdge(mapped_u, mapped_v, ew);
        });
    });

    return Gnew;
}

template <typename UnaryIdMapper>
Graph getRemappedGraph(const Graph& graph, count numNodes, UnaryIdMapper&& oldIdToNew, bool preallocate = true) {
    return getRemappedGraph(graph, numNodes,
        std::forward<UnaryIdMapper>(oldIdToNew), [](node) { return false; }, preallocate);
}


}	// namespace GraphTools
}	// namespace NetworKit

#endif // GRAPHTOOLS_H
