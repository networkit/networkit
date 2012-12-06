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

cluster Clustering::getCluster(node u) {
	return (*this)[u];
}

void Clustering::addToCluster(cluster c, node u) {
	assert ((*this)[u] == this->defaultValue);
	(*this)[u] = c;
}

void Clustering::toSingleton(node u) {
	(*this)[u] = this->nextCluster;
	this->nextCluster++;
}

void Clustering::moveToCluster(cluster c, node u) {
	assert ((*this)[u] != this->defaultValue);
	(*this)[u] = c;
}

void Clustering::mergeClusters(cluster c, cluster d) {
	cluster e = this->nextCluster;
	this->nextCluster++;
	for (node u = 1; u <= this->n; ++u) {
		if (( (*this)[u] == c) || (*this)[u] == d) {
			this->moveToCluster(e, u);
		}
	}
}

bool Clustering::isProper(const Graph& G) {
	// TODO: Clustering::isProper
}

} /* namespace EnsembleClustering */
