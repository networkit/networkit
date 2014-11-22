#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H

#include <unordered_map>
#include "../graph/Graph.h"

namespace NetworKit {

namespace GraphTools {

Graph getCompactedGraph(const Graph& graph);

std::unordered_map<node,node> getContinuousNodeIds(const Graph& graph);

Graph toUndirected(const Graph& graph);

}	// namespace GraphTools

}	// namespace NetworKit

#endif // GRAPHTOOLS_H