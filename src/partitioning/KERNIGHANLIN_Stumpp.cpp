/*
 * Louvain.cpp
 *
 *  Created on: 25.02.2013
 *      Author: Matthias Stumpp
 */

#include "KERNIGHANLIN_Stumpp.h"

namespace NetworKit {

cluster CL0 = 0;
cluster CL1 = 1;


KERNIGHAN_LIN::KERNIGHAN_LIN() {
}

KERNIGHAN_LIN::~KERNIGHAN_LIN() {
}

// compute initial gain for all nodes in graph
void KERNIGHAN_LIN::computeGainAllNodes(Graph& g, Clustering& c, NodeMap<double>& gain) {
	g.forNodes([&](node u) {
		this->computeGainForNode(g, u, c, gain);
	});
}

// compute initial gain for a certain node
void KERNIGHAN_LIN::computeGainForNode(Graph& g, node u, Clustering& c, NodeMap<double>& gain) {
	double samePart = 0.0, otherPart = 0.0;
	g.forWeightedNeighborsOf(u, [&](node v, edgeweight w) {
		if (c.inSameCluster(u, v)) {
			samePart += w;
		} else {
			otherPart += w;
		}
	});
	gain[u] = otherPart-samePart;
}

Clustering KERNIGHAN_LIN::partition(Graph& g, Clustering& c) {

	EdgeCut edgeCut;
	double cut = edgeCut.getQuality(c, g);
	DEBUG("KL start with edge cut " << cut);

	// unmark all nodes
	// false -> not visited, true otherwise
	NodeMap<char> marker(g.numberOfNodes(), 'f');
	// keeps the gain
	NodeMap<double> gain(g.numberOfNodes(), 0.0);
	// compute initial gain for all nodes
	TRACE("KL: start computing initial gains");
	this->computeGainAllNodes(g, c, gain);
	TRACE("KL: have computed initial gains");

	bool improved = false;
	std::vector<count> numCl = c.clusterSizes();
	count min = *min_element(std::begin(numCl), std::end(numCl));
	for (count n=0; n<min; n++) {
		TRACE("KL iteration " << n);

		double currentBest = 0.0;
		node best_u = none;
		node best_v = none;
		auto computeBestGain = [&](node u, node v, edgeweight w) {
			// continue if nodes in same cluster 
			// or u or v have been vistited already
			if (c.inSameCluster(u, v) || marker[u] == 't' || marker[v] == 't')
				return;

			// compute edge gain
			double runningGain = gain[u]+gain[v]-2*w;
			if (runningGain > currentBest) {
				improved = true;
				currentBest = runningGain;
				best_u = u;
				best_v = v;
			}
		};
		// FIXME: all node pairs from different blocks
		g.forWeightedEdges(computeBestGain);

		// no more improvements
		if (currentBest == 0.0)
			break;

		// otherwise, move nodes to other cluster
		if (c.clusterOf(best_u) == CL0) {
			c.moveToCluster(CL1, best_u);
			c.moveToCluster(CL0, best_v);
		} else {
			c.moveToCluster(CL0, best_u);
			c.moveToCluster(CL1, best_v);
		}
		// set marker to true, meaning node has been visited
		marker[best_u] = 't';
		marker[best_v] = 't';
		// recompute gain of nodes in neighbourhood of best two nodes
		// FIXME: iterate over neighbors of best_u and best_v
		this->computeGainForNode(g, best_u, c, gain);
		this->computeGainForNode(g, best_v, c, gain);
	}

	// either recursibely call or return
	if (improved)
		return this->partition(g, c);
	return c;
}

Clustering KERNIGHAN_LIN::run(Graph& graph) {
	Graph g = graph;

	// initialization
	// initial random partitioning
	ClusteringGenerator clusteringGenerator;
	TRACE("KL: start computing initial partition");
	Clustering clustering = clusteringGenerator.makeRandomClustering(g, 2);
	TRACE("KL: have computed initial partition, start KL pass");
	return this->partition(g, clustering);
}

std::string KERNIGHAN_LIN::toString() const {
	return "KERNIGHAN_LIN implementation by M. Stumpp";
}

} /* namespace NetworKit */
