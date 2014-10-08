/*
 *
 */

#include "UMST.h"

NetworKit::UMST::UMST(const NetworKit::Graph &G) : G(G) {

}

NetworKit::Graph NetworKit::UMST::generate() {
	Graph result(G.copyNodes());

	std::vector<weightedEdge<edgeweight> > weightedEdges;

	weightedEdges.reserve(G.numberOfEdges());

	G.forEdges([&](node u, node v, edgeweight weight) {
		weightedEdges.emplace_back(u, v, weight);
	});

	generate(std::move(weightedEdges), [&](weightedEdge<edgeweight> &e) {
		result.addEdge(e.u, e.v, e.attribute);
	});

	return result;
}
