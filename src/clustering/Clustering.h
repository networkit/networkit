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

	/**
	 * Construct new clustering.
	 *
	 * @param[in]	n	number of nodes
	 */
	Clustering(int64_t n);

	virtual ~Clustering();

	/**
	 *  Index operator.
	 *
	 *  @param[in]	u	a node
	 */
	inline cluster& operator [](const node& u) {
		return this->array[u];
	}
	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	inline const cluster& operator [](const node& u) const {
		return this->array[u];
	}

	/**
	 * Return the cluster (id) in which a node
	 * is contained.
	 */
	inline cluster& clusterOf(node u) {
		return (*this)[u];
	}

	/**
	 * Add a (previously unassigned) node to a cluster
	 */
	void addToCluster(cluster c, node u);

	/**
	 * Move a (previously assigned) node to a cluster.
	 */
	void moveToCluster(cluster c, node u);

	/**
	 * Creates a singleton cluster containing the node.
	 */
	void toSingleton(node u);

	/**
	 * Assigns the nodes from both clusters to a new cluster.
	 */
	void mergeClusters(cluster c, cluster d);

	/**
	 * Check whether this clustering is a proper clustering of
	 * the graph, i.e. a disjoint partition of the whole node set.
	 *
	 */
	bool isProper(const Graph& G);

	/**
	 * Get the lowest cluster id;
	 */
	cluster firstCluster();

	/**
	 * Get the highest cluster id that has been assigned.
	 * This gives an upper bound for the number of clusters in this clustering,
	 * although not the actual number of clusters since clusters can become empty.
	 */
	cluster lastCluster();

	/**
	 * Get iterator for all nodes in a cluster.
	 */
	// TODO: virtual Clustering::iterator iterCluster(cluster c);
};

} /* namespace EnsembleClustering */

#endif /* CLUSTERING_H_ */
