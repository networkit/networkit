/*
 * ClusteringGenerator2.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef ClusteringGenerator2_H_
#define ClusteringGenerator2_H_

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * Provides several methods for generating special clusterings.
 */
class ClusteringGenerator2 {

public:

	ClusteringGenerator2();

	virtual ~ClusteringGenerator2();

	/**
	 * Make a singleton clustering of G, i.e. a clustering in which every node
	 * belongs to its own cluster.
	 */
	virtual Partition makeSingletonClustering(Graph& G);

	/**
	 * Make a 1-clustering of G, i.e. a clustering in which all nodes belong to the same
	 * cluster.
	 */
	virtual Partition makeOneClustering(Graph& G);


	/**
	 * Make a clustering with k clusters to which the nodes are randomly assigned.
	 */
	virtual Partition makeRandomClustering(Graph& G, count k);


	/**
	 * Make a clustering with k clusters. The first n/k nodes are assigned to the
	 * first cluster, the next n/k nodes to the second cluster and so on.
	 */
	virtual Partition makeContinuousBalancedClustering(Graph& G, count k);

};

} /* namespace NetworKit */
#endif /* ClusteringGenerator2_H_ */
