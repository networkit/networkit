/*
 * Clustering.cpp
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#include "Clustering.h"

namespace EnsembleClustering {

Clustering::Clustering(int64_t n) : NodeMap<cluster>(n, 0), name("noname") {
	// all entries are initialized to 0, which means that the nodes are unclustered
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
	bool success = true;
	G.forallNodes([&](node v) {
		bool contained = this->contains(v);
		if (!contained) {
			success = false;
		}
	});
	return success;
}


int64_t Clustering::numberOfClusters() {
	int64_t k = 0; // number of clusters
	std::set<cluster> activeClusters;
	for (int64_t i = 1; i <= this->n; ++i) {
		cluster c = this->data[i];
		if (activeClusters.find(c) == activeClusters.end()) {
			k++;
			activeClusters.insert(c);
		}
	}
	return k;
}


cluster Clustering::upperBound() const {
	return this->n; // TODO: possibly inefficient
}


cluster Clustering::lowerBound() const {
	return 1;
}

void Clustering::allToSingletons() {
	#pragma omp parallel for
	for (node u = 1; u <= this->n; ++u) {
		this->data[u] = u;
	}
	this->nextCluster = this->n + 1;
}

cluster Clustering::addCluster() {
	cluster c = this->nextCluster;
	this->nextCluster += 1;
	return c;
}


bool Clustering::isInRange(node v) {
	return (1 <= v <= this->size());
}

bool Clustering::contains(node v) {
	assert (this->isInRange(v));	// assume that node is in range
	return (this->data[v] != 0);	// check if node is assigned to cluster, i.e. entry is not null
}

void Clustering::setName(std::string name) {
	this->name = name;
}

std::string Clustering::getName() const {
	return this->name;
}

void Clustering::print() const {
	std::cout << "{";
	for (int64_t i = 0; i <= this->n; ++i) {
		std::cout << i << ":" << this->data[i] << ", ";
	}
	std::cout << "}" << std::endl;
}

} /* namespace EnsembleClustering */


