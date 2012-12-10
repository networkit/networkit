/*
 * ClusteringGenerator.h
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#ifndef CLUSTERINGGENERATOR_H_
#define CLUSTERINGGENERATOR_H_

#include "Clustering.h"

namespace EnsembleClustering {

class ClusteringGenerator {

public:

	ClusteringGenerator();

	virtual ~ClusteringGenerator();

	/**
	 * Make a singleton clustering of G, i.e. a clustering in which every node
	 * belongs to its own cluster.
	 */
	virtual Clustering& makeSingletonClustering(const Graph& G);

	/**
	 * Make a 1-clustering of G, i.e. a clustering in which all nodes belong to the same
	 * cluster.
	 */
	virtual Clustering& makeOneClustering(const Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* CLUSTERINGGENERATOR_H_ */
