/*
 * DynamicBarabasiAlbertGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "../auxiliary/Random.h"

#include "DynamicBarabasiAlbertGenerator.h"

namespace NetworKit {


DynamicBarabasiAlbertGenerator::DynamicBarabasiAlbertGenerator(count k) : DynamicGraphSource(), k(k), degSum(0) {
	if (k <= 0) {
		throw std::runtime_error("k must be at least 1");
	}
}

void DynamicBarabasiAlbertGenerator::initializeGraph() {
	if (! this->graphSet) {
		throw std::runtime_error("Graph instance has not been set - call newGraph first");
	}

	// The network begins with an initial network of m0 nodes. m0 2 and the degree of each node in the initial network should be at least 1,
	// otherwise it will always remain disconnected from the rest of the network.
	for (count i = 0; i < k; ++i) {
		node u = Gproxy->addNode(); // assume node ids are assigned consecutively
		if (i > 0) {
			Gproxy->addEdge(u, u - 1); // connect to previous node to create a path
		}
	}

	degSum = 2 * G->numberOfEdges();

	this->graphInitialized = true; // graph has been properly initialized

	assert (G->numberOfNodes() == k);
	assert (G->numberOfEdges() == (k - 1));
}


void DynamicBarabasiAlbertGenerator::generate() {
	if (! this->graphInitialized) {
		throw std::runtime_error("Graph instance has not been initialized - call initializeGraph first");
	}

	// 3) go through the items one at a time, subtracting their weight from your random number, until you get the item where the random number is less than that item's weight
	node u = this->Gproxy->addNode();
	std::set<node> targets; // avoid duplicate edges by collecting target nodes in a set

	count nAttempts = 0;
	while (targets.size() < k) {
		nAttempts++;

//
//		if (nAttempts > 100) {
//			throw std::runtime_error("picking target nodes takes too long - something is wrong");
//		}

		// 2) pick a random number that is 0 or greater and is less than the sum of the weights
		uint64_t rand = Aux::Random::integer(degSum);

		bool found = false; // break from node iteration when done
		auto notFound = [&](){ return ! found; };
		this->G->forNodesWhile(notFound, [&](node v){
			if (v != u) { // skip u, which has degree 0 anyway, to avoid self-loops
				assert (rand >= 0);
				if (rand <= this->G->degree(v)) {
					found = true; // found a node to connect to
					targets.insert(v);
				}
				rand -= this->G->degree(v);
			}
		});
	}

	for (node v : targets) {
		this->Gproxy->addEdge(u, v);
		degSum += 2; 	// increment degree sum
	}

	this->Gproxy->timeStep(); // trigger a time step
}

} /* namespace NetworKit */
