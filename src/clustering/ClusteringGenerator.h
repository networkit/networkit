/*
 * ClusteringGenerator.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGGENERATOR_H_
#define CLUSTERINGGENERATOR_H_

#include "Clustering.h"

namespace NetworKit {

class ClusteringGenerator {

public:

	ClusteringGenerator();

	virtual ~ClusteringGenerator();

	/**
	 * Make a singleton clustering of G, i.e. a clustering in which every node
	 * belongs to its own cluster.
	 */
	virtual Clustering makeSingletonClustering(Graph& G);

	/**
	 * Make a 1-clustering of G, i.e. a clustering in which all nodes belong to the same
	 * cluster.
	 */
	virtual Clustering makeOneClustering(Graph& G);


	/**
	 * Make a clustering with k clusters to which the nodes are randomly assigned.
	 */
	virtual Clustering makeRandomClustering(Graph& G, count k);


	/**
	 * Make a clustering with k clusters. The first n/k nodes are assigned to the
	 * first cluster, the next n/k nodes to the second cluster and so on.
	 */
	virtual Clustering makeContinuousBalancedClustering(Graph& G, count k);

};

} /* namespace NetworKit */
#endif /* CLUSTERINGGENERATOR_H_ */
