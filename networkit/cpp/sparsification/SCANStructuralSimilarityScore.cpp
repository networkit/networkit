#include "SCANStructuralSimilarityScore.h"

NetworKit::SCANStructuralSimilarityScore::SCANStructuralSimilarityScore(const NetworKit::Graph &G, const std::vector< NetworKit::count > &triangles) : EdgeScore<double>(G), triangles(triangles) { }

void NetworKit::SCANStructuralSimilarityScore::run() {
	std::vector<double> scores(G.upperEdgeIdBound());

	if (!G.hasEdgeIds()) throw std::runtime_error("Error, edges must be indexed");

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		scores[eid] = (triangles[eid] + 1) * 1.0 / std::sqrt((G.degree(u) + 1)*(G.degree(v) + 1));
	});

	scoreData = std::move(scores);
	hasRun = true;
}
