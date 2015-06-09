/*
 * MatchingHierarchyBuilder.h
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MATCHINGHIERARCHYBUILDER_H_
#define MATCHINGHIERARCHYBUILDER_H_

#include "MultigridHierarchyBuilder.h"
#include "../algebraic/AdjacencyMatrix.h"
#include "../graph/Graph.h"
#include "../matching/PathGrowingMatcher.h"
#include "../coarsening/MatchingContracter.h"

namespace NetworKit {

class MatchingHierarchyBuilder : public MultigridHierarchyBuilder<MultigridHierarchy> {
private:
	const static count numCoarseGraphNodes = 30;

	template<typename InterpretationMatrix>
	MultigridHierarchy computeHierarchy(const Graph &graph) const;
	void interpolationMatrixFromMatching(const std::vector<node> &fineToCoarseIds, Matrix &interpolationMatrix) const;

public:
	MatchingHierarchyBuilder() {}

	MultigridHierarchy buildHierarchy(const Matrix &matrix) const;
	MultigridHierarchy buildHierarchy(const Graph &graph) const;
};

} /* namespace NetworKit */

#endif /* MATCHINGHIERARCHYBUILDER_H_ */
