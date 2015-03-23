/*
 * RandomNodeEdgeAttributizer.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include "RandomNodeEdgeAttributizer.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

RandomNodeEdgeAttributizer::RandomNodeEdgeAttributizer(const Graph& graph, double rneRatio) : graph(graph), rneRatio(rneRatio) {
}

std::vector< double > RandomNodeEdgeAttributizer::getAttribute() {
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	Graph backbone = graph;

	std::vector<double> edgeAttribute(graph.upperEdgeIdBound());

	count numRemoved = 0;

	std::vector< std::pair<node, node> > uniformlyRandomEdges;

	while (backbone.numberOfEdges() > 0) {
		if (Aux::Random::real() >= rneRatio) { // uniformly random
			bool edgeFound = false;

			while (!edgeFound) {
				if (uniformlyRandomEdges.empty()) {
					uniformlyRandomEdges = backbone.randomEdges(backbone.numberOfEdges() * (1.0 - rneRatio) + 20);
				}

				auto edge = uniformlyRandomEdges.back();
				uniformlyRandomEdges.pop_back();

				if (backbone.hasEdge(edge.first, edge.second)) {
					edgeid id = backbone.edgeId(edge.first, edge.second);

					edgeAttribute[id] = numRemoved * 1.0 / graph.numberOfEdges();

					backbone.removeEdge(edge.first, edge.second);

					edgeFound = true;
					++numRemoved;
				}
			}
		} else { // random node - edge
			auto edge = backbone.randomEdge();

			edgeid id = backbone.edgeId(edge.first, edge.second);

			edgeAttribute[id] = numRemoved * 1.0 / graph.numberOfEdges();

			backbone.removeEdge(edge.first, edge.second);

			++numRemoved;
		}
	}

	return edgeAttribute;
}

} /* namespace NetworKit */
