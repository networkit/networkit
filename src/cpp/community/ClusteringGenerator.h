/*
 * ClusteringGenerator.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGGENERATOR_H_
#define CLUSTERINGGENERATOR_H_

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * Provides several methods for generating special clusterings.
 */
class ClusteringGenerator {

public:

	/** Default destructor */
	virtual ~ClusteringGenerator();

	/**
	 * Make a singleton clustering of Graph @a G, i.e. a clustering in which every node
	 * belongs to its own cluster.
	 *
	 * @param G The graph.
	 * @return A Partition in which every node belongs to its own cluster.
	 */
	virtual Partition makeSingletonClustering(Graph& G);

	/**
	 * Make a 1-clustering of Graph @a G, i.e. a clustering in which all nodes belong to the same
	 * cluster.
	 *
	 * @param G The graph.
	 * @return A Partition in which all nodes belong to the same cluster.
	 */
	virtual Partition makeOneClustering(Graph& G);


	/**
	 * Make a clustering of Graph @a G with @a k clusters to which the nodes are randomly assigned.
	 *
	 * @param G The graph.
	 * @param k The amount of clusters.
	 * @return A Partition with @a k clusters and each node randomly assigned to one of them.
	 */
	virtual Partition makeRandomClustering(Graph& G, count k);


	/**
	 * Make a clustering of Graph @a G with @a k clusters. The first n/k nodes are assigned to the
	 * first cluster, the next n/k nodes to the second cluster and so on.
	 *
	 * @param G The graph.
	 * @param k The amount of clusters.
	 * @return A Partition with @a k clusters and each node assigned like described above.
	 */
	virtual Partition makeContinuousBalancedClustering(Graph& G, count k);

};

} /* namespace NetworKit */
#endif /* CLUSTERINGGENERATOR_H_ */
