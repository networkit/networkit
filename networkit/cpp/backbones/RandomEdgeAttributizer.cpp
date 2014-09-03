/*
 *
 */

#include "RandomEdgeAttributizer.h"
#include "../auxiliary/Random.h"

NetworKit::RandomEdgeAttributizer::RandomEdgeAttributizer(double rneRatio) : rneRatio(rneRatio) {

}

std::vector< double > NetworKit::RandomEdgeAttributizer::getAttribute(const NetworKit::Graph &g, const std::vector< int > &attribute) {
	Graph backbone = g;

	std::vector<double> edgeAttribute(g.upperEdgeIdBound());

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

					edgeAttribute[id] = numRemoved * 1.0 / g.numberOfEdges();

					backbone.removeEdge(edge.first, edge.second);

					edgeFound = true;
					++numRemoved;
				}
			}
		} else { // random node - edge
			auto edge = backbone.randomEdge();

			edgeid id = backbone.edgeId(edge.first, edge.second);

			edgeAttribute[id] = numRemoved * 1.0 / g.numberOfEdges();

			backbone.removeEdge(edge.first, edge.second);

			++numRemoved;
		}
	}

	return edgeAttribute;
}

