/*
 * MatchingHierarchyBuilder.cpp
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MatchingHierarchyBuilder.h"

namespace NetworKit {

template<typename InterpretationMatrix>
MultigridHierarchy MatchingHierarchyBuilder::computeHierarchy(const Graph &graph) const {
	Matrix mat = InterpretationMatrix(graph);
	MultigridHierarchy hierarchy(mat);
	PathGrowingMatcher matcher;
	MatchingContracter contracter;

	Graph currentGraph = graph;
	INFO("Computing hierarchy");
	do {
		// calculate matching of graph
		Matching matching = matcher.run(currentGraph, false);

		std::pair<Graph,std::vector<node>> coarsening = contracter.run(currentGraph, matching, true);

		// compute interpolation matrix
		Matrix interpolationMatrix(currentGraph.upperNodeIdBound(), coarsening.first.upperNodeIdBound());
		interpolationMatrixFromMatching(coarsening.second, interpolationMatrix);

		// add level to hierarchy
		currentGraph = coarsening.first;
		hierarchy.addLevel(InterpretationMatrix(currentGraph), interpolationMatrix);
		INFO("add level with nC=", currentGraph.numberOfNodes(), " and matched nodes ", matching.size());
	} while (currentGraph.upperNodeIdBound() > numCoarseGraphNodes);
	INFO("Done with ", hierarchy.getNumLevels(), " levels");

	return hierarchy;
}

void MatchingHierarchyBuilder::interpolationMatrixFromMatching(const std::vector<node> &fineToCoarseIds, Matrix &interpolationMatrix) const {
	for (index i = 0; i < fineToCoarseIds.size(); ++i) {
		interpolationMatrix.setValue(i, fineToCoarseIds[i], 1);
	}
}

MultigridHierarchy MatchingHierarchyBuilder::buildHierarchy(const Matrix &matrix) const {
	Graph graph(std::max(matrix.numberOfRows(), matrix.numberOfColumns()), true, true);
	matrix.forNonZeroElementsInRowOrder([&](node i, node j, edgeweight w){
		graph.addEdge(i, j, w);
	});

	return computeHierarchy<AdjacencyMatrix>(graph);
}

MultigridHierarchy MatchingHierarchyBuilder::buildHierarchy(const Graph &graph) const {
	return computeHierarchy<LaplacianMatrix>(graph);
}

} /* namespace NetworKit */
