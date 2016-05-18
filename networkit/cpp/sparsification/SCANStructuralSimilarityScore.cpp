#include "SCANStructuralSimilarityScore.h"

NetworKit::SCANStructuralSimilarityScore::SCANStructuralSimilarityScore(const NetworKit::Graph &G, const std::vector< NetworKit::count > &triangles) : EdgeScore<double>(G), triangles(triangles) { }

void NetworKit::SCANStructuralSimilarityScore::run() {
	std::vector<double> workScores(G.upperEdgeIdBound());

	if (!G.hasEdgeIds()) throw std::runtime_error("Error, edges must be indexed");

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		workScores[eid] = (triangles[eid] + 1) * 1.0 / std::sqrt((G.degree(u) + 1)*(G.degree(v) + 1));
	});

	scoreData = std::move(workScores);
	hasRun = true;
}

double NetworKit::SCANStructuralSimilarityScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double NetworKit::SCANStructuralSimilarityScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}
