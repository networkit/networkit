#include "../graph/Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup community
 */
namespace GraphClusteringTools {

float getImbalance(Partition& zeta);

Graph communicationGraph(const Graph& graph, Partition& zeta);

count weightedDegreeWithCluster(const Graph& graph, Partition& zeta, node u, index cid);

bool isProperClustering(Graph& G, Partition& zeta);

bool isSingletonClustering(Graph& G, Partition& zeta);

bool isOneClustering(Graph& G, Partition& zeta);

bool equalClusterings(Partition& zeta, Partition& eta, Graph& G);

}

}


