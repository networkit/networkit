/*
* EffectiveDiameter.cpp
*
*  Created on: 16.06.2014
*      Author: Marc Nemes
*/

#include "EffectiveDiameter.h"
#include "ConnectedComponents.h"
#include "../auxiliary/Random.h"

#include <math.h>
#include <iterator>
#include <stdlib.h>
#include <omp.h>
#include <map>

namespace NetworKit {

double EffectiveDiameter::effectiveDiameter(const Graph& G, const double ratio, const count k, const count r) {
	// the length of the bitmask where the number of connected nodes is saved
	count lengthOfBitmask = (count) ceil(log2(G.numberOfNodes()));
	// saves all k bitmasks for every node of the current iteration
	std::vector<std::vector<unsigned int> > mCurr;
	// saves all k bitmasks for every node of the previous iteration
	std::vector<std::vector<unsigned int> > mPrev;
	// the list of nodes that are already connected to all other nodes
	std::vector<node> finishedNodes;
	// the maximum possible bitmask based on the random initialization of all k bitmasks
	std::vector<count> highestCount;
	// the amount of nodes that need to be connected to all others nodes
	count threshold = (count) (ceil(ratio * G.numberOfNodes()));
	// the current distance of the neighborhoods
	count h = 1;
	// sums over the number of edges needed to reach 90% of all other nodes
	double effectiveDiameter = 0;
	// the estimated number of connected nodes
	double estimatedConnectedNodes;
	// used for setting a random bit in the bitmasks
	double random;
	// the position of the bit that has been set in a bitmask
	count position;
	// nodes that are not connected to enough nodes yet
	std::vector<node> activeNodes;

	// initialize all vectors
	highestCount.assign(k, 0);
	G.forNodes([&](node v) {
		finishedNodes.push_back(v);
		finishedNodes[v] = 0;
		std::vector<unsigned int> bitmasks;
		bitmasks.assign(k, 0);
		mCurr.push_back(bitmasks);
		mPrev.push_back(bitmasks);
		activeNodes.push_back(v);
		// set one bit in each bitmask with probability P(bit i=1) = 0.5^(i+1), i=0,..
		for (count j = 0; j < k; j++) {
			random = Aux::Random::real(0,1);
			position = ceil(log(random)/log(0.5) - 1);
			// set the bit in the bitmask
			if (position < lengthOfBitmask+r) {
				mPrev[v][j] |= 1 << position;
			}
			// add the current bit to the maximum-bitmask
			highestCount[j] = highestCount[j] | mPrev[v][j];
		}
	});

	// as long as we need to connect more nodes
	while (!activeNodes.empty()) {
		for (count x = 0; x < activeNodes.size(); x++) {
			node v = activeNodes[x];
			#pragma omp parallel for
			// for each parallel approximation
			for (count j = 0; j < k; j++) {
				// the node is still connected to all previous neighbors
				mCurr[v][j] = mPrev[v][j];
				// and to all previous neighbors of all its neighbors
				G.forNeighborsOf(v, [&](node u) {
					mCurr[v][j] = mCurr[v][j] | mPrev[u][j];
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

			// calculate the estimated number of neighbors where 0.77351 is a correction factor and the result of a complex sum
			estimatedConnectedNodes = (pow(2,b) / 0.77351);

			// check whether all k bitmask for this node have reached their highest possible value
			bool nodeFinished = true;
			for (count j = 0; j < k; j++) {
				if (mCurr[v][j] != highestCount[j]) {
					nodeFinished = false;
					break;
				}
			}
			// if the node wont change or is connected to enough nodes it must no longer be considered
			if (estimatedConnectedNodes >= threshold || nodeFinished) {
				effectiveDiameter += h;
				// remove the current node from future iterations
				std::swap(activeNodes[x], activeNodes.back());
				activeNodes.pop_back();
			}
		}
		mPrev = mCurr;
		h++;
	}
	return effectiveDiameter/G.numberOfNodes();
}

double EffectiveDiameter::effectiveDiameterExact(const Graph& G, const double ratio) {
	// saves the reachable nodes of the current iteration
	std::vector<std::vector<bool> > mCurr;
	// saves the reachable nodes of the previous iteration
	std::vector<std::vector<bool> > mPrev;
	// sums over the number of edges needed to reach 90% of all other nodes
	double effectiveDiameter = 0;
	// the current distance of the neighborhoods
	count h = 1;
	// number of nodes that need to be connected with all other nodes
	count threshold = (uint64_t) (ceil(ratio * G.numberOfNodes()) + 0.5);
	// nodes that are not connected to enough nodes yet
	std::vector<node> activeNodes;

	// initialize all nodes
	G.forNodes([&](node v){
		std::vector<bool> connectedNodes;
		// initialize n entries with value 0
		connectedNodes.assign(G.upperNodeIdBound(),0);
		// the node is always connected to itself
		connectedNodes[v] = 1;
		mCurr.push_back(connectedNodes);
		mPrev.push_back(connectedNodes);
		activeNodes.push_back(v);
	});

	// as long as we need to connect more nodes
	while (!activeNodes.empty()) {
		for (count x = 0; x < activeNodes.size(); x++) {
			node v = activeNodes[x];
				mCurr[v] = mPrev[v];
				G.forNeighborsOf(v, [&](node u) {
					for (count i = 0; i < G.numberOfNodes(); i++) {
						// add the current neighbor of u to the neighborhood of v
						mCurr[v][i] = mCurr[v][i] || mPrev[u][i];
					}
				});

				// compute the number of connected nodes
				count numConnectedNodes = 0;
				for (count i = 0; i < G.numberOfNodes(); i++) {
					if (mCurr[v][i] == 1) {
						numConnectedNodes++;
					}
				}

				// when the number of connected nodes surpasses the threshold the node must no longer be considered
				if (numConnectedNodes >= threshold) {
					effectiveDiameter += h;
					// remove the current node from future iterations
					std::swap(activeNodes[x], activeNodes.back());
					activeNodes.pop_back();
					x--; //don't skip former activeNodes.back() that has been switched to activeNodes[x]
				}
			}
			mPrev = mCurr;
			h++;
		}
	// return the found effective diameter
	return effectiveDiameter/G.numberOfNodes();
}

std::map<count, double> EffectiveDiameter::hopPlot(const Graph& G, const count maxDistance, const count k, const count r) {
	//the returned hop-plot
	std::map<count, double> hopPlot;
	// the length of the bitmask where the number of connected nodes is saved
	count lengthOfBitmask = (count) ceil(log2(G.numberOfNodes()));
	// saves all k bitmasks for every node of the current iteration
	std::vector<std::vector<unsigned int> > mCurr;
	// saves all k bitmasks for every node of the previous iteration
	std::vector<std::vector<unsigned int> > mPrev;
	// the maximum possible bitmask based on the random initialization of all k bitmasks
	std::vector<count> highestCount;
	// the current distance of the neighborhoods
	count h = 1;
	// the estimated number of connected nodes
	double estimatedConnectedNodes;
	// the sum over all estimated connected nodes
	double totalConnectedNodes;
	// used for setting a random bit in the bitmasks
	double random;
	// the position of the bit that has been set in a bitmask
	count position;
	// nodes that are not connected to enough nodes yet
	std::vector<node> activeNodes;

	// initialize all vectors
	highestCount.assign(k, 0);
	G.forNodes([&](node v) {
		std::vector<unsigned int> bitmasks;
		bitmasks.assign(k, 0);
		mCurr.push_back(bitmasks);
		mPrev.push_back(bitmasks);
		activeNodes.push_back(v);

		// set one bit in each bitmask with probability P(bit i=1) = 0.5^(i+1), i=0,..
		for (count j = 0; j < k; j++) {
			random = Aux::Random::real(0,1);
			position = ceil(log(random)/log(0.5) - 1);
			// set the bit in the bitmask
			if (position < lengthOfBitmask+r) {
				mPrev[v][j] |= 1 << position;
			}
			// add the current bit to the maximum-bitmask
			highestCount[j] = highestCount[j] | mPrev[v][j];
		}
	});
	// at zero distance, all nodes can only reach themselves
	hopPlot[0] = 1/G.numberOfNodes();
	// as long as we need to connect more nodes
	while (!activeNodes.empty() && (maxDistance <= 0 || h < maxDistance)) {
		totalConnectedNodes = 0;
		for (count x = 0; x < activeNodes.size(); x++) {
			node v = activeNodes[x];
			// for each parallel approximation
			for (count j = 0; j < k; j++) {
				// the node is still connected to all previous neighbors
				mCurr[v][j] = mPrev[v][j];
				// and to all previous neighbors of all its neighbors
				G.forNeighborsOf(v, [&](node u) {
					mCurr[v][j] = mCurr[v][j] | mPrev[u][j];
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
			// For the origin of the factor 0.77351 see http://www.mathcs.emory.edu/~cheung/papers/StreamDB/Probab/1985-Flajolet-Probabilistic-counting.pdf Theorem 3.A (p. 193)
			estimatedConnectedNodes = (pow(2,b) / 0.77351);

			// enforce monotonicity
			if (estimatedConnectedNodes > G.numberOfNodes()) {
				estimatedConnectedNodes = G.numberOfNodes();
			}

			// check whether all k bitmask for this node have reached the highest possible value
			bool nodeFinished = true;
			for (count j = 0; j < k; j++) {
				if (mCurr[v][j] != highestCount[j]) {
					nodeFinished = false;
					break;
				}
			}

			// if the node wont change or is connected to enough nodes it must no longer be considered
			if (estimatedConnectedNodes >= G.numberOfNodes() || nodeFinished) {
				// remove the current node from future iterations
				std::swap(activeNodes[x], activeNodes.back());
				activeNodes.pop_back();
				totalConnectedNodes += G.numberOfNodes();
			} else {
				// add value of the node to all nodes so we can calculate the average
				totalConnectedNodes += estimatedConnectedNodes;
			}
		}
		// add nodes that are already connected to all nodes
		totalConnectedNodes += (G.numberOfNodes() - activeNodes.size()) * G.numberOfNodes();
		// compute the fraction of connected nodes
		hopPlot[h] = totalConnectedNodes/(G.numberOfNodes()*G.numberOfNodes());
		if (hopPlot[h] > 1) {
			hopPlot[h] = 1;
		}
		mPrev = mCurr;
		h++;
	}
	return hopPlot;
}

}
