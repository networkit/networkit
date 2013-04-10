/*
 * DynamicBarabasiAlbertGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "DynamicBarabasiAlbertGenerator.h"

namespace NetworKit {

DynamicBarabasiAlbertGenerator::DynamicBarabasiAlbertGenerator() {
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


	this->initializeGraph();

	count degSum = 2; // TODO: parameters


	while (! terminate()) {

		// 1) calculate the sum of all the weights
		// TODO: update degSum


		// 2) pick a random number that is 0 or greater and is less than the sum of the weights

		int k = 2; // number of edges per new node

		int64_t rand; // TODO: generate random number from [0, degSum]

		// 3) go through the items one at a time, subtracting their weight from your random number, until you get the item where the random number is less than that item's weight
		node u = this->proxy->addNode();

		for (int i = 0; i < k; ++i) {
			bool found = false; // break from node iteration when done
			auto notFound = [&](){ return ! found; };
			this->G->forNodesWhile(notFound, [&](node v){
				if (rand < this->G->degree(v)) {
					found = true; // found a node to connect to
					this->proxy->addEdge(u, v);
					degSum += 2; 	// increment degree sum
				}
			});
		}

	}





}

} /* namespace NetworKit */
