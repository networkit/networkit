/*
 * MultiLevelSetup.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef MULTILEVELSETUP_H_
#define MULTILEVELSETUP_H_

#include "../../algebraic/Matrix.h"
#include "LevelHierarchy.h"
#include "../Smoother.h"

#include "../../algebraic/CSRMatrix.h"

namespace NetworKit {

#define UNDECIDED -1

class MultiLevelSetup {

private:
	const Smoother &smoother;

	bool coarseningElimination(CSRMatrix &matrix, LevelHierarchy &hierarchy) const;
	count lowDegreeSweep(const CSRMatrix &matrix, std::vector<bool> &fNode, index stage) const;
	void eliminationOperators(const CSRMatrix &matrix, const std::vector<index> &fSet, const std::vector<index> &coarseIndex, CSRMatrix &P, Vector &q) const;
	void subMatrix(const CSRMatrix &matrix, const std::vector<index> &rows, const std::vector<index> &columns, const std::vector<index> &coarseIndex, CSRMatrix &result) const;

	void coarseningAggregation(CSRMatrix &matrix, LevelHierarchy &hierarchy, Vector &tv, count numTVVectors) const;
	std::vector<Vector> generateTVs(const CSRMatrix &matrix, Vector &tv, const count numVectors) const;
	void addHighDegreeSeedNodes(const CSRMatrix &matrix, std::vector<int64_t> &status) const;
	void aggregateLooseNodes(const CSRMatrix &strongAdjMatrix, std::vector<int64_t> &status, count &nc) const;
	void computeStrongAdjacencyMatrix(const CSRMatrix &matrix, CSRMatrix &strongAdjMatrix) const;
	void computeAffinityMatrix(const CSRMatrix &matrix, const std::vector<Vector> &tVs, CSRMatrix &affinityMatrix) const;
	void aggregationStage(const CSRMatrix &matrix, count &nc, const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, std::vector<Vector> &tVs, std::vector<int64_t> &status) const;
	void computeStrongNeighbors(const CSRMatrix &affinityMatrix, const std::vector<int64_t> &status, std::vector<std::vector<index>> &bins) const;
	bool findBestSeed(const CSRMatrix &affinityMatrix, const std::vector<index> &strongNeighborsOfU, const std::vector<int64_t> &status, const index u, index &s) const;
	bool findBestSeedEnergyCorrected(const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, const std::vector<double> &diag, const std::vector<Vector> &tVs, const std::vector<int64_t> &status, const index u, index &s) const;

	bool canCoarsen(const CSRMatrix &A) const;
	bool isRelaxationFast(const CSRMatrix &A, index lvlIndex, Vector &tv) const;

	void galerkinOperator(const CSRMatrix &P, const CSRMatrix &A, const std::vector<index> &PColIndex, const std::vector<std::vector<index>> &PRowIndex, CSRMatrix &B) const;

public:
	MultiLevelSetup(const Smoother &smoother);

	void setup(const Graph &G, LevelHierarchy &hierarchy) const;
	void setup(const CSRMatrix &matrix, LevelHierarchy &hierarchy) const;
};

} /* namespace NetworKit */

#endif /* MULTILEVELSETUP_H_ */
