/*
 * Clustering.h
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <set>

#include "../../graph/NodeMap.h"


namespace EnsembleClustering {

typedef int64_t cluster;	//!< cluster is represented as a 1-based index

class Clustering: public NodeMap<cluster> {

protected:

	cluster nextCluster;	//!< next free cluster id for new cluster

	inline cluster getNextCluster() {
		// TODO: performance - is this a bottleneck?
		cluster c = this->nextCluster;
		this->nextCluster++;
		return c;
	}

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
		return this->data[u];
	}
	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	inline const cluster& operator [](const node& u) const {
		return this->data[u];
	}

	/**
	 * Return the cluster (id) in which a node
	 * is contained.
	 */
	inline cluster clusterOf(node u) const {
		return this->data[u];
	}

	/**
	 * Call this before assigning nodes to cluster ids.
	 */
	cluster addCluster();

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
	 * Assigns every node to a singleton cluster.
	 */
	void allToSingletons();

	/**
	 * Assigns the nodes from both clusters to a new cluster.
	 */
	void mergeClusters(cluster c, cluster d);

	/**
	 * Check whether this clustering is a proper clustering of
	 * the graph, i.e. a disjoint partition of the whole node set.
	 *
	 */
	bool isProper(Graph& G);


	/**
	 * Get the current number of clusters in this clustering.
	 */
	int64_t numberOfClusters();


	/**
	 * Return an upper bound for the cluster ids that have been assigned.
	 */
	cluster upperBound() const;

	/**
	 * Return a lower bound for the cluster ids that have been assigned.
	 */
	cluster lowerBound() const;


	// DEBUG

	void print();


};

} /* namespace EnsembleClustering */

#endif /* CLUSTERING_H_ */
