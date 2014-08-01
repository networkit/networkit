#include "../graph/Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup community
 */
namespace GraphClusteringTools {

float getImbalance(const Partition& zeta);

Graph communicationGraph(const Graph& graph, Partition &zeta);

count weightedDegreeWithCluster(const Graph& graph, const Partition& zeta, node u, index cid);

bool isProperClustering(const Graph& G, const Partition& zeta);

bool isSingletonClustering(const Graph &G, const Partition &zeta);

bool isOneClustering(const Graph& G, const Partition& zeta);

bool equalClusterings(const Partition &zeta, const Partition &eta, Graph& G);

}

}


