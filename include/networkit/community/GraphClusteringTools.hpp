/*
 * GraphClusteringTools.hpp
 */

#ifndef NETWORKIT_COMMUNITY_GRAPH_CLUSTERING_TOOLS_HPP_
#define NETWORKIT_COMMUNITY_GRAPH_CLUSTERING_TOOLS_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup community
 */
namespace GraphClusteringTools {

/**
 * Get the imbalance of clusters in the given partition.
 *
 * @param zeta Input partition
 * @return total imbalance between clusters
 */
float getImbalance(const Partition &zeta);

/**
 * Get the imbalance of clusters in the given partition compared to a given graph.
 *
 * @param zeta Input partition
 * @param graph Input graph
 * @return total imbalance between clusters
 */
float getImbalance(const Partition &zeta, const Graph &graph);

/**
 * Get the communication graph for a given graph and its partition.
 * A communication graph consists of a number of nodes, which equal
 * the number of clusters in the partition. The edges between nodes
 * in the communication graph account for the total edge weight for all
 * edges between two clusters. For unweighted graphs, the edge weight in
 * the communication graph is equal to the number of edges between two
 * clusters.
 *
 * @param graph The input graph
 * @param zeta Partition, which contains information about clusters in the graph
 * @return communication graph
 */
Graph communicationGraph(const Graph &graph, Partition &zeta);

/**
 * Get weightedDegree of node u for a cluster (represented by a partition) of index cid.
 *
 * @param graph The input graph
 * @param zeta Partition, which contains information about clusters in the graph
 * @param u node
 * @param cid index of cluster
 * @return weighted degree of node u for cluster index cid
 */
count weightedDegreeWithCluster(const Graph &graph, const Partition &zeta, node u, index cid);

/**
 * Check whether a partition is a proper clustering for a given graph
 *
 * @param graph The input graph
 * @param zeta Partition, which contains information about clusters in the graph
 * @return True if the partition is a proper clustering, False if not
 */
bool isProperClustering(const Graph &G, const Partition &zeta);

/**
 * Check whether a partition is a proper singleton clustering for a given graph
 *
 * @param graph The input graph
 * @param zeta Partition, which contains information about clusters in the graph
 * @return True if the partition is a proper singleton clustering, False if not
 */
bool isSingletonClustering(const Graph &G, const Partition &zeta);

/**
 * Check whether a partition is a proper one clustering for a given graph
 *
 * @param graph The input graph
 * @param zeta Partition, which contains information about clusters in the graph
 * @return True if the partition is a proper one clustering, False if not
 */
bool isOneClustering(const Graph &G, const Partition &zeta);

/**
 * Check whether two paritions are equal for a given graph
 *
 * @param graph The input graph
 * @param zeta Partition, which contains information about clusters in the graph
 * @param eta Partition, which contains information about clusters in the graph
 * @return True if both partitions are the same, False if not
 */
bool equalClusterings(const Partition &zeta, const Partition &eta, Graph &G);

} // namespace GraphClusteringTools

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_GRAPH_CLUSTERING_TOOLS_HPP_
