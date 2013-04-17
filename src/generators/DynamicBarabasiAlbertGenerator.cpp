/*
 * DynamicBarabasiAlbertGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "DynamicBarabasiAlbertGenerator.h"

namespace NetworKit {

DynamicBarabasiAlbertGenerator::DynamicBarabasiAlbertGenerator(GraphEventProxy& proxy, count k) : DynamicGraphGenerator(proxy), k(k) {
	// TODO Auto-generated constructor stub

}

void DynamicBarabasiAlbertGenerator::initializeGraph() {
	// The network begins with an initial network of m0 nodes. m0 ³ 2 and the degree of each node in the initial network should be at least 1,
	// otherwise it will always remain disconnected from the rest of the network.
	for (count i = 0; i < k; ++i) {
		node u = Gproxy->addNode(); // assume first node has id 0
		if (u > 0) {
			Gproxy->addEdge(u, u - 1); // connect to previous node to create a path
		}
	}
	assert (G->numberOfNodes() == k);
	assert (G->numberOfEdges() == (k - 1));
}

DynamicBarabasiAlbertGenerator::~DynamicBarabasiAlbertGenerator() {
	// TODO Auto-generated destructor stub
}


void DynamicBarabasiAlbertGenerator::generateWhile(std::function<bool(void)> cont) {
	INFO("[BEGIN] generating graph");

	assert (G->numberOfNodes() >= k); // there must be at least as many nodes in the graphs as the number of edges added in each step
	count degSum = 2; // TODO: parameters
	Aux::RandomInteger randInt;

	while (cont()) {

		// 3) go through the items one at a time, subtracting their weight from your random number, until you get the item where the random number is less than that item's weight
		node u = this->Gproxy->addNode();
		std::set<node> targets; // avoid duplicate edges by collecting target nodes in a set

		while (targets.size() < k) {
			TRACE("pick random node to connect to");
			DEBUG("targets.size() == " << targets.size());
			// 2) pick a random number that is 0 or greater and is less than the sum of the weights
			int64_t rand = randInt.generate(0, degSum);

			bool found = false; // break from node iteration when done
			auto notFound = [&](){ return ! found; };
			this->G->forNodesWhile(notFound, [&](node v){
				if (v != u) { // skip u, which has degree 0 anyway, to avoid self-loops
					assert (rand >= 0);
					if (rand < this->G->degree(v)) {
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

	} // end while

	INFO("[STOP]Êgenerating graph");

}

} /* namespace NetworKit */
