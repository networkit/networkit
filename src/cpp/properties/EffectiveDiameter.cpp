/*
 * EffectiveDiameter.cpp
 *
 *  Created on: 16.06.2014
 *      Author: Marc Nemes
 */

#include "EffectiveDiameter.h"
#include "ConnectedComponents.h"

#include <math.h>
#include <iterator>
#include <stdlib.h>

namespace NetworKit {

count EffectiveDiameter::effectiveDiameter(const Graph& G) {
	return EffectiveDiameter::effectiveDiameter(G, 0.9, 64, 7, 4);
}

count EffectiveDiameter::effectiveDiameter(const Graph& G, const double ratio) {
	return EffectiveDiameter::effectiveDiameter(G, ratio, 64, 7, 4);
}

/*
this is a variaton of the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining
in Massive Graphs" by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
*/
count EffectiveDiameter::effectiveDiameter(const Graph& G, const double ratio, const count k, const count r, const count l) {
	// check whether the graph is connected
	ConnectedComponents cc(G);
	cc.run();
	if (cc.numberOfComponents() > 1) {
		throw std::runtime_error("Graph not connected - diameter is infinite");
	}

	// the length of the bitmask where the number of connected nodes is saved
	count lengthOfBitmask = (count) ceil(log2(G.numberOfNodes()));
	// saves all k bitmasks for every node of the current iteration
	std::vector<std::vector<unsigned int> > mCurr;
	// saved all k bitmasks for every node of the previous iteration
	std::vector<std::vector<unsigned int> > mLast;
	// the list of nodes that are already connected to all other nodes
	std::vector<node> finishedNodes;
	// the amount of estimated connected nodes from the previous iteration
	std::vector<double> previousEstimatedConnectedNodes;
	// the amount of iterations that a node passed without connecting to more nodes
	std::vector<count> consecutiveCounter;
	// the amount of nodes that need to be connected to all others nodes
	count threshold = (count) (ceil(ratio * G.numberOfNodes()));
	// the current distance of the neighborhoods
	count h = 1;
	// the amount of nodes that are connected to all other nodes (|finishedNodes|)
	count numberOfFinishedNodes = 0;
	double random;
	srand (time(NULL));

	// initialize all vectors
	G.forNodes([&](node v) {
		finishedNodes.push_back(v);
		finishedNodes[v] = 0;
		previousEstimatedConnectedNodes.push_back(v);
		previousEstimatedConnectedNodes[v] = -1;
		consecutiveCounter.push_back(v);
		consecutiveCounter[v] = 0;
		std::vector<unsigned int> tmp;
		for (count j = 0; j < k; j++) {
			tmp.push_back(0);
		}
		mCurr.push_back(tmp);
		mLast.push_back(tmp);

		// set one bit in each bitmask with probability P(bit i=1) = 0.5^(i+1), i=0,..
		for (count j = 0; j < k; j++) {
			random = (rand()/(double)(RAND_MAX));
			for (count i = 0; i < lengthOfBitmask+r; i++) {
				if (random > pow(0.5,i+1)) {
					mLast[v][j] |= 1 << i;
					break;
				}
			}
		}
	});

	// as long as we need to connect more nodes
	while (numberOfFinishedNodes < threshold) {
		G.forNodes([&](node v) {
			// if the current node is not yet connected to all other nodes
			if (finishedNodes[v] == 0) {
				// for each parallel approximation
				for (count j = 0; j < k; j++) {
					// the node is still connected to all previous neighbors
					mCurr[v][j] = mLast[v][j];
					// and to all previous neighbors of all its neighbors
					G.forNeighborsOf(v, [&](node u) {
						mCurr[v][j] = mCurr[v][j] | mLast[u][j];
					});
				}
				// the estimated number of neighbors
				double estimatedConnectedNodes;
				// the least bit number in the bitmask of the current node/distance that has not been set
				double b = 0;

				for (count j = 0; j < k; j++) {
					for (count i = 0; i < sizeof(i)*8; i++) {
						if (((mCurr[v][j] >> i) & 1) == 0) {
							b += i;
							break;
						}
					}
				}

				// calculate the average least bit number that has not been set over all parallel approximations
				b = b / k;
				// calculate the estimated number of neighbors
				estimatedConnectedNodes = (pow(2,b) / 0.77351);
				// when the amount of neighbors is still the same increase the counter so we now when this node must no longer be considered
				if (previousEstimatedConnectedNodes[v] == estimatedConnectedNodes) {
					consecutiveCounter[v]++;
				} else {
					consecutiveCounter[v] = 0;
				}
				// this node is probably connected to all other nodes and must no longer be considered
				if (consecutiveCounter[v] >= l) {
					finishedNodes[v] = 1;
					numberOfFinishedNodes++;
				}
				previousEstimatedConnectedNodes[v] = estimatedConnectedNodes;
			}
		});
		mLast = mCurr;
		h++;
	}
	// we always mark a node l iterations too late
	h = h-l-1;
	return h;
}

count EffectiveDiameter::effectiveDiameterExact(const Graph& G) {
	// Check whether the graph is connected
	ConnectedComponents cc(G);
	cc.run();
	if (cc.numberOfComponents() > 1) {
		throw std::runtime_error("Graph not connected - diameter is infinite");
	}

	// list of nodes that are connected to all other nodes
	std::set<node> finished;
	// diameter[node][distance][connected_nodes]
	std::vector<std::vector<std::set<node> > > diameter;
	// initialize all nodes
	G.forNodes([&](node v){
		std::set<node> connectedNodes;
		std::vector<std::set<node> > inner;
		// at the beginning, nodes are only connected to themselves
		connectedNodes.insert(v);
		// connect all nodes with themselves at the distance 0
		inner.push_back(connectedNodes);
		diameter.push_back(inner);
	});

	// current diameter
	count d = 1;
	// number of nodes that are connected to all other nodes
	count i = 0;
	// number of nodes that need to be connected with all other nodes
	count threshold = (uint64_t) (ceil(0.9 * G.numberOfNodes()) + 0.5);
	// as long as we need to connect more nodes
	while (i < threshold) {
		G.forNodes([&](node v){
			// only consider nodes that are not already connected to every other node
			if (finished.find(v) == finished.end()) {
				// the node is connected to all nodes from the previous iteration
				std::set<node> connectedNodes = diameter[v][d-1];
				diameter[v].push_back(connectedNodes);
				// and to all previous connected nodes of all neighbors
				G.forNeighborsOf(v, [&](node u) {
					for (auto it : diameter[u][d-1]) {
						// add the current neighbor of u to the neighborhood of v
						diameter[v][d].insert(it);
					}
				});
				// do no longer consider this node once it's connected to all nodes
				if (diameter[v][d].size() == G.numberOfNodes()) {
					finished.insert(v);
					i++;
				}
			}
		});
		d++;
	}
	// return the found effective diameter
	return d-1;
}

}
