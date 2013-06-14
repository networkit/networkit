/*
 * DynamicBarabasiAlbertGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "DynamicBarabasiAlbertGenerator.h"

namespace NetworKit {

DynamicBarabasiAlbertGenerator::DynamicBarabasiAlbertGenerator() {
}

DynamicBarabasiAlbertGenerator::DynamicBarabasiAlbertGenerator(GraphEventProxy& proxy, count k) : DynamicGraphGenerator(proxy), k(k) {
	if (k <= 0) {
		throw std::runtime_error("k must be at least 1");
	}
	if (this->G->numberOfNodes() > 0) {
		throw std::runtime_error("this generator needs to be initialized with an empty graph");
	}
}

void DynamicBarabasiAlbertGenerator::initializeGraph() {
	// The network begins with an initial network of m0 nodes. m0 ³ 2 and the degree of each node in the initial network should be at least 1,
	// otherwise it will always remain disconnected from the rest of the network.
	DEBUG("k = " << k);
	for (int i = 0; i < k; ++i) {
		TRACE("adding node");
		node u = Gproxy->addNode(); // assume node ids are assigned consecutively
		if (i > 0) {
			Gproxy->addEdge(u, u - 1); // connect to previous node to create a path
		}
	}

	degSum = 2 * G->numberOfEdges();

	DEBUG("n = " << G->numberOfNodes());
	assert (G->numberOfNodes() == k);
	assert (G->numberOfEdges() == (k - 1));
}

DynamicBarabasiAlbertGenerator::~DynamicBarabasiAlbertGenerator() {
	// TODO Auto-generated destructor stub
}


void DynamicBarabasiAlbertGenerator::generate() {
	// 3) go through the items one at a time, subtracting their weight from your random number, until you get the item where the random number is less than that item's weight
	node u = this->Gproxy->addNode();
	std::set<node> targets; // avoid duplicate edges by collecting target nodes in a set

	count nAttempts = 0;
	while (targets.size() < k) {
		nAttempts++;

		TRACE("pick random node to connect to");
		TRACE("targets.size() == " << targets.size());

		if (nAttempts > 100) {
			throw std::runtime_error("picking target nodes takes too long - something is wrong");
		}

		// 2) pick a random number that is 0 or greater and is less than the sum of the weights
		int64_t rand = randInt.generate(0, degSum);

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
