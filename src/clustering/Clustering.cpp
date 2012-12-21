/*
 * Clustering.cpp
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#include "Clustering.h"

namespace EnsembleClustering {

Clustering::Clustering(int64_t n) : NodeMap<cluster>(n, 0) {
	this->nextCluster = 1; //!< first cluster index is 1
}

Clustering::~Clustering() {
	// TODO Auto-generated destructor stub
}



void Clustering::addToCluster(cluster c, node u) {
	assert((*this)[u] == this->defaultValue);
	(*this)[u] = c;
}

void Clustering::toSingleton(node u) {
	(*this)[u] = this->getNextCluster();
}

void Clustering::moveToCluster(cluster c, node u) {
	assert((*this)[u] != this->defaultValue);
	(*this)[u] = c;
}

void Clustering::mergeClusters(cluster c, cluster d) {
	cluster e = this->getNextCluster();
	for (node u = 1; u <= this->n; ++u) {
		if (((*this)[u] == c) || (*this)[u] == d) {
			this->moveToCluster(e, u);
		}
	}
}



bool Clustering::isProper(Graph& G) {
	// test whether each node has been assigned to a cluster
	G.forallNodes([&](node v) {
		cluster c = this->clusterOf(v);
		assert (c <= this->upperBound());
	});
	bool success = true;
	return success;
}


int64_t Clustering::numberOfClusters() {
	int64_t k = 0; // number of clusters
	std::set<cluster> activeClusters;
	for (int64_t i = 1; i <= this->n; ++i) {
		cluster c = this->array[i];
		if (activeClusters.find(c) == activeClusters.end()) {
			k++;
			activeClusters.insert(c);
		}
	}
	return k;
}


cluster Clustering::upperBound() const {
	return this->n;
}

void Clustering::allToSingletons() {
	#pragma omp parallel for
	for (node u = 1; u <= this->n; ++u) {
		this->array[u] = u;
	}
	this->nextCluster = this->n + 1;
}

cluster Clustering::addCluster() {
	cluster c = this->nextCluster;
	this->nextCluster += 1;
	return c;
}

cluster Clustering::lowerBound() const {
	return 1;
}

void Clustering::print() {
	std::cout << "{";
	for (int64_t i = 0; i <= this->n; ++i) {
		std::cout << i << ":" << this->array[i] << ", ";
	}
	std::cout << "}" << std::endl;
}

} /* namespace EnsembleClustering */


