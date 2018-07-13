/*
 * MocnikGeneratorBasic.cpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#include "MocnikGeneratorBasic.h"
#include "../auxiliary/Random.h"
#include <unordered_map>

namespace NetworKit {

MocnikGeneratorBasic::MocnikGeneratorBasic(count dim, count n, double k): dim(dim), n(n), k(k) {
}

// GEOMETRY

// norm of a vector
static inline double norm(std::vector<double> &v, const double &shift) {
	double x = 0;
	for (count j = 0; j < v.size(); j++) {
		x += (v[j] + shift) * (v[j] + shift);
	}
	return sqrt(x);
}

// Euclidean distance between two vectors
static inline double dist(std::vector<double> &v, std::vector<double> &w) {
	double x = 0;
	for (count j = 0; j < v.size(); j++) {
		x += std::pow(v[j] - w[j], 2);
	}
	x = sqrt(x);
	return x;
}

// GRAPH GENERATION

Graph MocnikGeneratorBasic::generate() {
	assert (dim > 0);
	assert (n > 0);
	assert (k > 1);

	// vector containing the nodes resp. their positions
	NodePositionMap nodesPosition;

	// create graph
	Graph G(0, false, true);

	// create the nodes
	node curr = 0;
	while (curr < n) {
		std::vector<double> v = {};
		for (count j = 0; j < dim; j++) {
			v.push_back(Aux::Random::real());
		}
		// test wheather the new node would be contained in the ball B_{.5}(.5, ..., .5)
		if (norm(v, -.5) < .5) {
			nodesPosition[G.addNode()] = v;
			curr++;
		}
	}

	// create the edges
	double x;
	for (auto kv : nodesPosition) {
		// compute the minimal distance from a node to all other nodes
		double dist_min = -1;
		for (auto kv2 : nodesPosition) {
			x = dist(kv.second, kv2.second);
			if (kv.first != kv2.first && (x < dist_min || dist_min == -1)) {
				dist_min = x;
			}
		}
		// add the edges
		for (auto kv2 : nodesPosition) {
			if (dist(kv.second, kv2.second) <= k * dist_min && kv.first != kv2.first) {
				G.addEdge(kv.first, kv2.first);
			}
		}
	}

	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
