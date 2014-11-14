/*
 * LAMG.h
 *
 *  Created on: 12.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LAMG_H_
#define LAMG_H_

#include "MultigridHierarchyBuilder.h"
#include "LAMGHierarchy.h"
#include "../algebraic/Matrix.h"
#include "../algebraic/DiagonalMatrix.h"
#include "../graph/Graph.h"
#include "../auxiliary/Random.h"
#include "GaussSeidelRelaxation.h"

namespace NetworKit {

class LAMG : public MultigridHierarchyBuilder {
private:
	GaussSeidelRelaxation smoother;

	std::vector<bool> lowDegreeNodes(const Matrix &matrix);
	std::vector<Vector> computeTestVectors(const Matrix &matrix, const count numberOfVectors) const;
	Matrix computeAffinityMatrix(const Matrix &matrix, const std::vector<Vector> &testVectors) const;

	bool addEliminationLevel(Matrix &matrix, LAMGHierarchy &hierarchy);

	bool addAggregationLevel(Matrix &matrix, LAMGHierarchy &hierarchy);

public:
	LAMG();

	LAMGHierarchy buildHierarchy(const Matrix &matrix);
	LAMGHierarchy buildHierarchy(const Graph &graph);

};

} /* namespace NetworKit */

#endif /* LAMG_H_ */
