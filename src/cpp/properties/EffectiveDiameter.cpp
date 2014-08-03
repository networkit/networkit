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
#include <omp.h>

namespace NetworKit {

double EffectiveDiameter::effectiveDiameter(const Graph& G) {
	return EffectiveDiameter::effectiveDiameter(G, 0.9, 64, 7);
}

double EffectiveDiameter::effectiveDiameter(const Graph& G, const double ratio) {
	return EffectiveDiameter::effectiveDiameter(G, ratio, 64, 7);
}

/*
this is a variaton of the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining
in Massive Graphs" by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
*/
double EffectiveDiameter::effectiveDiameter(const Graph& G, const double ratio, const count k, const count r) {
	// set number of threads for OpenMP
	omp_set_num_threads(4);
	// the length of the bitmask where the number of connected nodes is saved
	count lengthOfBitmask = (count) ceil(log2(G.numberOfNodes()));
	// saves all k bitmasks for every node of the current iteration
	std::vector<std::vector<unsigned int> > mCurr;
	// saves all k bitmasks for every node of the previous iteration
	std::vector<std::vector<unsigned int> > mLast;
	// the list of nodes that are already connected to all other nodes
	std::vector<node> finishedNodes;
	// the maximum possible bitmask based on the random initialization of all k bitmasks
	std::vector<count> highestCount;
	// the amount of nodes that need to be connected to all others nodes
	count threshold = (count) (ceil(ratio * G.numberOfNodes()));
	// the current distance of the neighborhoods
	count h = 1;
	// the amount of nodes that are connected to all other nodes (|finishedNodes|)
	count numberOfFinishedNodes = 0;
	// the minimal number of leading ones in all bitmasks
	count leadingOnes = 0; // TODO
	// sums over the number of edges needed to reach 90% of all other nodes
	double effectiveDiameter = 0;
	// the estimated number of connected nodes
	double estimatedConnectedNodes;

	double random;
	srand (time(NULL));

	// initialize all vectors
	for (count j = 0; j < k; j++) {
		highestCount.push_back(j);
		highestCount[j] = 0;
	}
	G.forNodes([&](node v) {
		finishedNodes.push_back(v);
		finishedNodes[v] = 0;
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
			// add the current bit to the maximum-bitmask
			highestCount[j] = highestCount[j] | mLast[v][j];
		}
	});

	// as long as we need to connect more nodes
	while (numberOfFinishedNodes < G.numberOfNodes()) {
		G.forNodes([&](node v) {
			// if the current node is not yet connected to all other nodes
			if (finishedNodes[v] == 0) {
				#pragma omp parallel for
				// for each parallel approximation
				for (count j = 0; j < k; j++) {
					// the node is still connected to all previous neighbors
					mCurr[v][j] = mLast[v][j];
					// and to all previous neighbors of all its neighbors
					G.forNeighborsOf(v, [&](node u) {
						mCurr[v][j] = mCurr[v][j] | mLast[u][j];
					});
				}
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

				bool nodeFinished = true;
				for (count j = 0; j < k; j++) {
					if (mCurr[v][j] != highestCount[j]) {
						nodeFinished = false;
						break;
					}
				}
				if (estimatedConnectedNodes >= threshold || nodeFinished) {
					finishedNodes[v] = 1;
					numberOfFinishedNodes++;
					effectiveDiameter += h;
				}
			}
		});
		mLast = mCurr;
		h++;
	}
	return effectiveDiameter/G.numberOfNodes();
}

count EffectiveDiameter::effectiveDiameterExact(const Graph& G) {
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
