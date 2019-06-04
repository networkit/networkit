#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H

#include <unordered_map>
#include <networkit/graph/Graph.hpp>

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
 * @param numNodes Number of nodes in the output graph.
 * @param oldIdToNew Translate old id to new ones.
 */
template <typename UnaryFunc>
Graph remapNodes(const Graph& graph, count numNodes, UnaryFunc oldIdToNew) {
    // TODO: Add reserveNeighbors as soon as Graph supports it
#ifndef NDEBUG
    graph.forNodes([&] (node u) {assert(oldIdToNew(u) < numNodes);});
#endif

    Graph Gnew(numNodes, graph.isWeighted(), graph.isDirected());

    graph.forEdges([&](const node u, const node v, edgeweight ew) {
        Gnew.addEdge(oldIdToNew(u), oldIdToNew(v), ew);
    });

    return Gnew;
}


}	// namespace GraphTools

}	// namespace NetworKit

#endif // GRAPHTOOLS_H

