/*
 * RmatGenerator.cpp
 *
 *  Created on: 18.03.2014
 *      Author: Henning, cls
 */

#include "RmatGenerator.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/NumericTools.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

RmatGenerator::RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d, bool weighted, count reduceNodes):
	scale(scale), edgeFactor(edgeFactor), a(a), b(b), c(c), d(d), weighted(weighted), reduceNodes(reduceNodes)
{
    if (scale > 63) throw std::runtime_error("Cannot generate more than 2^63 nodes");
	double sum = a+b+c+d;
	INFO("sum of probabilities: ", sum);
	if (!Aux::NumericTools::equal(sum, 1.0, 0.0001)) throw std::runtime_error("Probabilities in Rmat have to sum to 1.");
	defaultEdgeWeight = 1.0;
}

Graph RmatGenerator::generate() {
	count n = (1 << scale);
	count numEdges = n * edgeFactor;
	Graph G(n, true);
	double ab = a+b;
	double abc = ab+c;

	auto quadrant([&]() {
		double r = Aux::Random::probability();
		TRACE("r: ", r);

		if (r <= a) {
			return 0;
		}
		else if (r <= ab) {
			return 1;
		}
		else if (r <= abc) {
			return 2;
		}
		else return 3;
	});

	auto drawEdge([&]() {
		node u = 0;
		node v = 0;
		for (index i = 0; i < scale; ++i) {
			count q = quadrant();
//			TRACE("q: ", q);
			u = u << 1;
			v = v << 1;
			u = u | (q >> 1);
			v = v | (q & 1);
		}

		return std::make_pair(u, v);
	});

	for (index e = 0; e < numEdges; ++e) {
		std::pair<node, node> drawnEdge = drawEdge();
//		TRACE("edge drawn: ", drawnEdge.first, " - ", drawnEdge.second);
		G.increaseWeight(drawnEdge.first, drawnEdge.second, defaultEdgeWeight);
	}

	// delete random nodes to achieve node count
	INFO("deleting random nodes: ", reduceNodes);
	for (count i = 0; i < reduceNodes; ++i) {
		node u = G.randomNode();
		std::vector<std::pair<node, node>> incidentEdges;
		G.forEdgesOf(u, [&](node u, node v) {
			incidentEdges.push_back({u,v});
		});
		for (auto edge : incidentEdges) {
			node x, y;
			std::tie(x, y) = edge;
			G.removeEdge(x, y);
		}
		assert (G.degree(u) == 0);
		G.removeNode(u);
	}

	if (!weighted) {
		// set unit weights
		G.forEdges([&](node u, node v) {
			G.setWeight(u, v, 1.0);
		});
	}

	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
