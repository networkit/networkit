/*
 * Clustering.h
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include "../graph/NodeMap.h"

namespace EnsembleClustering {

typedef int64_t cluster;	//!< cluster is represented as a 1-based index

class Clustering: public NodeMap<cluster> {

protected:

	cluster nextCluster;	//!< next free cluster id for new cluster

public:

	Clustering(int64_t n);

	virtual ~Clustering();

	cluster getCluster(node u);

	void addToCluster(cluster c, node u);

	void toSingleton(node u);

	void moveToCluster(cluster c, node u);

	void mergeClusters(cluster c, cluster d);

	/**
	 * Check whether this clustering is a proper clustering of
	 * the graph, i.e. a disjoint partition of the whole node set.
	 *
	 */
	virtual bool isProper(const Graph& G);

};

} /* namespace EnsembleClustering */
#endif /* CLUSTERING_H_ */
