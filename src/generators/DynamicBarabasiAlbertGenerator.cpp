/*
 * DynamicBarabasiAlbertGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "DynamicBarabasiAlbertGenerator.h"

namespace NetworKit {

DynamicBarabasiAlbertGenerator::DynamicBarabasiAlbertGenerator(GraphEventProxy& proxy, int k) : DynamicGraphGenerator(proxy), k(k) {
	// TODO Auto-generated constructor stub

}

void DynamicBarabasiAlbertGenerator::initializeGraph() {
	// The network begins with an initial network of m0 nodes. m0 ³ 2 and the degree of each node in the initial network should be at least 1,
	// otherwise it will always remain disconnected from the rest of the network.
	// TODO: make n_0 a parameter, here it is fixed to 2
	node s1 = this->proxy->addNode();
	node s2 = this->proxy->addNode();
	this->proxy->addEdge(s1, s2);
}

DynamicBarabasiAlbertGenerator::~DynamicBarabasiAlbertGenerator() {
	// TODO Auto-generated destructor stub
}


void DynamicBarabasiAlbertGenerator::generate(std::function<bool(void)> terminate) {


	this->initializeGraph(); // TODO: not here when generation needs to be stopped and resumed

	count degSum = 2; // TODO: parameters

	Aux::RandomInteger randInt;


	while (! terminate()) {


		// 3) go through the items one at a time, subtracting their weight from your random number, until you get the item where the random number is less than that item's weight
		node u = this->proxy->addNode();
		TRACE("adding node " << u);

		for (int i = 0; i < this->k; ++i) {
			TRACE("i = " << i);
			// 2) pick a random number that is 0 or greater and is less than the sum of the weights
			int64_t rand = randInt.generate(0, degSum); // TODO: generate random number from [0, degSum]

			bool found = false; // break from node iteration when done
			auto notFound = [&](){ return ! found; };
			this->G->forNodesWhile(notFound, [&](node v){
				TRACE("scanning node " << v);
				rand -= this->G->degree(v);
				TRACE("rand = " << rand);
				if (rand < this->G->degree(v)) {
					found = true; // found a node to connect to
					this->proxy->addEdge(u, v);
					degSum += 2; 	// increment degree sum
				}
			});
		}

	} // end while

}

} /* namespace NetworKit */
