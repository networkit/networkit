/*
 * Lamg.h
 *
 *  Created on: Oct 20, 2015
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_NUMERICS_LAMG_LAMG_H_
#define NETWORKIT_CPP_NUMERICS_LAMG_LAMG_H_

#include <vector>

#include "../LinearSolver.h"
#include "MultiLevelSetup.h"
#include "SolverLamg.h"
#include "../GaussSeidelRelaxation.h"

namespace NetworKit {

class Lamg : public LinearSolver {
private:
	bool validSetup;
	GaussSeidelRelaxation smoother;
	MultiLevelSetup lamgSetup;
	CSRMatrix laplacianMatrix;
	std::vector<LevelHierarchy> compHierarchies;
	std::vector<SolverLamg> compSolvers;
	std::vector<LAMGSolverStatus> compStati;

	std::vector<Vector> initialVectors;
	std::vector<Vector> rhsVectors;

	count numComponents;
	std::vector<std::vector<index>> components;
	std::vector<index> graph2Components;

	void initializeForOneComponent();

public:
	Lamg(const double desiredResidualRed = 1e-6);
	~Lamg() = default;

	void setup(const CSRMatrix &laplacianMatrix);
	void setupConnected(const CSRMatrix &laplacianMatrix);

	SolverStatus solve(const Vector &rhs, Vector &result, count maxConvergenceTime = 5 * 60 * 1000, count maxIterations = std::numeric_limits<count>::max());

};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_NUMERICS_LAMG_LAMG_H_ */
