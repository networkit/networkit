#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H

#include <unordered_map>
#include "../graph/Graph.h"

namespace NetworKit {

namespace GraphTools {


/**
 * Computes a graph with the same structure but with continuous node ids.
 * @param  graph     The graph to be compacted.
 * @param  nodeIdMap The map providing the information about the node ids.
 * @return           Returns a compacted Graph.
 */
Graph getCompactedGraph(const Graph& graph, std::unordered_map<node,node>& nodeIdMap);

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
std::vector<node> invertContinuousNodeIds(std::unordered_map<node,node>& nodeIdMap, const Graph& G);

/**
 * Constructs a new graph that has the same node ids as before it was compacted.
 * @param  invertedIdMap The node id mapping from continuous node ids to noncontinuous node ids.
 * @param  G             The compacted graph.
 * @return               The original graph.
 */
Graph restoreGraph(std::vector<node>& invertedIdMap, const Graph& G);




}	// namespace GraphTools

}	// namespace NetworKit

#endif // GRAPHTOOLS_H