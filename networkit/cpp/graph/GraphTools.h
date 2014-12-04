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

}	// namespace GraphTools

}	// namespace NetworKit

#endif // GRAPHTOOLS_H