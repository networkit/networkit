/*
 * Clustering.cpp
 *
 *  Created on: 31.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Clustering.h"
#include <algorithm>

namespace NetworKit {

Clustering::Clustering() : NodeMap<cluster>(0, none), name("noname")  {
	this->nextCluster = 0; // first cluster index is 0
	this->upperIdBound = 0; // upper id bound = n is okay only for agglomeratively created clusters
}


Clustering::Clustering(count n) :
		NodeMap<cluster>(n, none), name("noname") {
	// all entries are initialized to none, which means that the nodes are unclustered
	this->nextCluster = 0; // first cluster index is 0
	this->upperIdBound = n; // upper id bound = n is okay only for agglomeratively created clusters
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
	for (node u = 0; u < this->n; ++u) {
		if (((*this)[u] == c) || (*this)[u] == d) {
			this->moveToCluster(e, u);
		}
	}
}

bool Clustering::isProper(Graph& G) {
	// test whether each node has been assigned to a cluster
	bool success = true;
	G.forNodes([&](node v) {
		bool contained = this->contains(v);
		if (!contained) {
			ERROR("Clustering does not contain node " << v);
			success = false;
		}
	});
	return success;
}

count Clustering::numberOfClusters() const {

	count nc = this->upperBound();
	std::vector<int> exists(nc, 0); // a boolean vector would not be thread-safe

	this->parallelForEntries([&](node u, cluster C) {
		if (C != none) {
			exists[C] = 1;
		}
	});

	count k = 0; // number of actually existing clusters
	#pragma omp parallel for reduction(+:k)
	for (index i = 0; i < nc; ++i) {
		if (exists[i]) {
			k++;
		}
	}

	return k;
}

cluster Clustering::upperBound() const {
	return this->upperIdBound;
}

cluster Clustering::lowerBound() const {
	return 0;
}

void Clustering::allToSingletons() {
	DEBUG("this->n is " << this->n);
	for (node u = 0; u < this->n; ++u) {
		this->data[u] = u;
	}
	this->nextCluster = this->n;
}

cluster Clustering::addCluster() {
	cluster c = this->nextCluster;
	this->nextCluster += 1;
	return c;
}

bool Clustering::isInRange(node v) {
	return ((0 <= v) && (v < this->numberOfNodes()));
}

bool Clustering::contains(node v) {
	// assert (this->isInRange(v));	// assume that node is in range
	return (this->data[v] != this->defaultValue); // check if node is assigned to cluster, i.e. entry is not null
}

void Clustering::setName(std::string name) {
	this->name = name;
}

std::string Clustering::getName() const {
	return this->name;
}

bool Clustering::isOneClustering(Graph& G) {
	return (this->numberOfClusters() == 1);
}

bool Clustering::isSingletonClustering(Graph& G) {
	return (this->numberOfClusters() == G.numberOfNodes());
}

bool Clustering::inSameCluster(node u, node v) {
	assert(this->contains(u));
	assert(this->contains(v));
	return (this->clusterOf(u) == this->clusterOf(v));
}

bool Clustering::equals(Clustering& other, Graph& G) {
	bool eq = true;
	G.parallelForEdges([&](node u, node v) {
		if (this->inSameCluster(u, v)) {
			if (!other.inSameCluster(u, v)) {
				eq = false;
			}
		}
		else {
			if (other.inSameCluster(u, v)) {
				eq = false;
			}
		}

	});
	return eq;
}

void Clustering::setUpperBound(cluster id) {
	this->upperIdBound = id;
}


void Clustering::compact() {
	std::map<cluster, cluster> compactingMap;
	cluster nextIndex = 0;

	this->forEntries([&](node v, cluster c) {
		if (compactingMap.count(c) == 0) {
			// insert and increase nextIndex
			if (c != none) {
				compactingMap.insert(std::make_pair(c, nextIndex));
				++nextIndex;
			}
			else {
				compactingMap.insert(std::make_pair(none, none));
			}
		}
	});

	this->parallelForEntries([&](node v, cluster c) {
		data[v] = compactingMap[c];
	});

	setUpperBound(nextIndex);
	TRACE("upperBound: " << upperBound());
}

std::vector<count> Clustering::clusterSizes() const{
	count numC = this->numberOfClusters();
	std::vector<count> clusterSizes(numC);

	this->parallelForEntries([&](node v, cluster c) {
		c = data[v];
		++clusterSizes[c];
	});


	return clusterSizes;
}


void Clustering::append(node u) {
	this->data.push_back(this->defaultValue);
	this->n += 1;
	assert (this->data[u] == this->defaultValue); // assumption: push_back creates entry at index u
}

} /* namespace NetworKit */

