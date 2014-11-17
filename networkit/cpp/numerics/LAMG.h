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

#include <limits>
#include <unordered_map>

namespace NetworKit {

class LAMG : public MultigridHierarchyBuilder {
private:
	GaussSeidelRelaxation smoother;
	double guard, cycleIndex;

	std::vector<bool> lowDegreeNodes(const Matrix &matrix) const;
	std::vector<Vector> computeTestVectors(const Matrix &matrix, const count numOfVectors) const;
	Matrix computeAffinityMatrix(const Matrix &matrix, const std::vector<Vector> &testVectors) const;
	void addHighDegreeSeedNodes(const Matrix &matrix, std::vector<int64_t> &status) const;
	std::vector<std::vector<index>> computeStrongNeighbors(const Matrix &affinityMatrix, const double delta) const;
	int64_t findBestSeed(const Matrix &affinityMatrix, const std::vector<index> &strongNeighbors, const std::vector<int64_t> &status, const index u) const;
	void aggregationStage(const Matrix &matrix, count &numCoarseNodes, const Matrix &affinityMatrix, std::vector<Vector> &testVectors, std::vector<int64_t> &status, std::vector<count> &aggregateSize, double delta) const;

	bool addEliminationLevel(Matrix &matrix, LAMGHierarchy &hierarchy);

	bool addAggregationLevel(Matrix &matrix, LAMGHierarchy &hierarchy, const count numTestVectors);

public:
	LAMG(double guard=0.7, double cycleIndex=1.5);

	LAMGHierarchy buildHierarchy(const Matrix &matrix);
	LAMGHierarchy buildHierarchy(const Graph &graph);

};

} /* namespace NetworKit */

#endif /* LAMG_H_ */
