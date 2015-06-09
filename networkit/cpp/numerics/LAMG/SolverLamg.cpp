/*
 * SolverLamg.cpp
 *
 *  Created on: 12.01.2015
 *      Author: Michael
 */

#include "SolverLamg.h"
#include "LAMGSettings.h"
//#include "Level.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include "../../io/LineFileReader.h"

#include "../../auxiliary/Enforce.h"
#include "../../auxiliary/Timer.h"

namespace NetworKit {


SolverLamg::SolverLamg(LevelHierarchy &hierarchy, const Smoother &smoother) : hierarchy(hierarchy), smoother(smoother) {
}

void SolverLamg::solve(Vector &x, const Vector &b, LAMGSolverStatus &status) {
	if (hierarchy.size() >= 2) {
		Vector bc = b;
		Vector xc = x;
		int finest = 0;

		if (hierarchy.getType(1) == ELIMINATION) {
			hierarchy.at(1).restrict(b, bc);
			if (hierarchy.at(1).getLaplacian().numberOfRows() == 1) {
				x = 0.0;
				return;
			} else {
				hierarchy.at(1).coarseType(x, xc);
				finest = 1;
			}
		}

		solveCycle(xc, bc, finest, status);

		if (finest == 1) { // interpolate from finest == ELIMINATION level back to actual finest level
			hierarchy.at(1).interpolate(xc, x);
		} else {
			x = xc;
		}
	} else {
		solveCycle(x, b, 0, status);
	}

	double residual = (b - hierarchy.at(0).getLaplacian() * x).length();
	status.residual = residual;
	DEBUG("final residual\t ", residual);
}

void SolverLamg::solveCycle(Vector &x, const Vector &b, int finest, LAMGSolverStatus &status) {
	Aux::Timer timer;
	timer.start();

	// data structures for iterate recombination
	history = std::vector<std::vector<Vector>>(hierarchy.size());
	rHistory = std::vector<std::vector<Vector>>(hierarchy.size());
	latestIterate = std::vector<index>(hierarchy.size(), 0);
	numActiveIterates = std::vector<count>(hierarchy.size(), 0);
	for (index i = 0; i < hierarchy.size(); ++i) {
		history[i] = std::vector<Vector>(MAX_COMBINED_ITERATES, Vector(hierarchy.at(i).getNumberOfNodes()));
		rHistory[i] = std::vector<Vector>(MAX_COMBINED_ITERATES, Vector(hierarchy.at(i).getNumberOfNodes()));
	}

	Vector r = b - hierarchy.at(finest).getLaplacian() * x;
	double residual = r.length();
	double finalResidual = residual * status.desiredResidual;
	double lastResidual = std::numeric_limits<double>::max();

	count iterations = 0;
	status.residualHistory.emplace_back(residual);
	while (residual > finalResidual && residual < lastResidual && iterations < status.maxIters && timer.elapsedMilliseconds() <= status.maxConvergenceTime ) {
		DEBUG("iter ", iterations, " r=", residual);
		lastResidual = residual;

		cycle(x, b, finest, status);
		r = b - hierarchy.at(finest).getLaplacian() * x;
		residual = r.length();
		status.residualHistory.emplace_back(residual);
		iterations++;
	}

	timer.stop();

	status.numIters = iterations;
	status.residual = r.length();
	status.converged = r.length() <= finalResidual;
	DEBUG("nIter\t ", iterations);
}

void SolverLamg::cycle(Vector &x, const Vector &b, int finest, const LAMGSolverStatus &status) {
	std::vector<count> numVisits(hierarchy.size()-1, 0);
	int coarsest = hierarchy.size() - 1;
	std::vector<Vector> X(hierarchy.size());
	std::vector<Vector> B(hierarchy.size());
	X[finest] = x;
	B[finest] = b;




	int currLvl = finest;
	int nextLvl = 0;
	double maxVisits = 0.0;

	saveIterate(currLvl, X[currLvl], B[currLvl] - hierarchy.at(currLvl).getLaplacian() * X[currLvl]);

	while (true) {
		if (currLvl == coarsest) {
			nextLvl = currLvl - 1;
			if (currLvl == finest) { // finest level
				X[currLvl] = smoother.relax(hierarchy.at(currLvl).getLaplacian(), B[currLvl], X[currLvl], 1);
			} else {
				X[currLvl] = smoother.relax(hierarchy.at(currLvl).getLaplacian(), B[currLvl], X[currLvl], SETUP_RELAX_COARSEST_SWEEPS);
			}
		} else {
			if (currLvl == finest) {
				maxVisits = 1.0;
			} else {
				maxVisits = hierarchy.cycleIndex(currLvl) * numVisits[currLvl-1];
			}

			if (numVisits[currLvl] < maxVisits) {
				nextLvl = currLvl + 1;
			} else {
				nextLvl = currLvl - 1;
			}
		}

		if (nextLvl < finest) break;

		if (nextLvl > currLvl) {  // preProcess
			numVisits[currLvl]++;

			Vector xf = X[currLvl];
			Vector bf = B[currLvl];
			Vector bc;
			Vector xc;
			if (hierarchy.getType(nextLvl) != ELIMINATION) {
				xf = smoother.relax(hierarchy.at(currLvl).getLaplacian(), bf, xf, status.numPreSmoothIters);
				X[currLvl] = xf;
			}

			if (hierarchy.getType(nextLvl) == ELIMINATION) {
				hierarchy.at(nextLvl).restrict(bf, bc);
			} else {
				Vector rf = bf - hierarchy.at(currLvl).getLaplacian() * xf;
				hierarchy.at(nextLvl).restrict(rf, bc);
			}

			hierarchy.at(nextLvl).coarseType(xf, xc);
			B[nextLvl] = bc;
			X[nextLvl] = xc;
			clearHistory(nextLvl);
			if (currLvl > finest) {
				clearHistory(currLvl);
			}
		} else { // postProcess
			Vector xf = X[nextLvl];
			Vector bf = B[nextLvl];
			Vector rf = bf - hierarchy.at(nextLvl).getLaplacian() * xf;
			Vector xc = X[currLvl];
			Vector rc = B[currLvl] - hierarchy.at(currLvl).getLaplacian() * xc;

			if (currLvl == coarsest || hierarchy.getType(currLvl+1) != ELIMINATION) {
				minRes(currLvl, xc, rc);
				X[currLvl] = xc;
			}

			if (nextLvl > finest) {
				saveIterate(nextLvl, xf, rf);
			}


			if (hierarchy.getType(currLvl) == ELIMINATION) {
				hierarchy.at(currLvl).interpolate(xc, xf);
			} else {
				hierarchy.at(currLvl).interpolate(xc, xf);
				xf += X[nextLvl];
			}

			if (hierarchy.getType(currLvl) != ELIMINATION) {
				xf = smoother.relax(hierarchy.at(nextLvl).getLaplacian(), bf, xf, status.numPostSmoothIters);
			}

			X[nextLvl] = xf;
		}

		currLvl = nextLvl;
	} // while


	// post-cycle finest
	if (hierarchy.size() > finest + 1 && hierarchy.getType(finest+1) != ELIMINATION) { // do an iterate recombination on calculated solutions
//		DEBUG("numActive: ", numActiveIterates[finest]);
		minRes(finest, X[finest], B[finest] - hierarchy.at(finest).getLaplacian() * X[finest]);
	}




	double mean = X[finest].mean();
#pragma omp parallel for
	for (index i = 0; i < X[finest].getDimension(); ++i) { // subtract mean from all entries in final solution vector
		X[finest][i] -= mean;
	}

	x = X[finest];
}

void SolverLamg::multigridCycle(index level, Vector &xf, const Vector &bf) {
	if (level == hierarchy.size()-1) { // at coarsest level
		if (level == 0) {
			xf = smoother.relax(hierarchy.at(level).getLaplacian(), bf, xf, 1);
		} else {
			xf = smoother.relax(hierarchy.at(level).getLaplacian(), bf, xf, SETUP_RELAX_COARSEST_SWEEPS);
		}

		return;
	}

	// pre-smoothing
	if (hierarchy.getType(level+1) != ELIMINATION) {
		xf = smoother.relax(hierarchy.at(level).getLaplacian(), bf, xf, 1);
	}

	// restriction
	Vector bc;
	Vector xc = xf;
	if (hierarchy.getType(level+1) == AGGREGATION) { // coarse b is residual
		hierarchy.at(level+1).restrict(bf - hierarchy.at(level).getLaplacian() * xf, bc);
	} else {
		hierarchy.at(level+1).restrict(bf, bc);
	}

	hierarchy.at(level+1).coarseType(xf, xc);

	// solve coarse problem
	multigridCycle(level+1, xc, bc);

	// interpolation
	Vector newXf = xf;
	hierarchy.at(level+1).interpolate(xc, newXf);
	if (hierarchy.getType(level+1) == AGGREGATION) {
		newXf += xf;
	}

	xf = newXf;

	// post-smoothing
	if (hierarchy.getType(level+1) != ELIMINATION) {
		xf = smoother.relax(hierarchy.at(level).getLaplacian(), bf, xf, 2);
	}
}

void SolverLamg::saveIterate(index level, const Vector &x, const Vector &r) {
	// update latest pointer
	index i = latestIterate[level];
	latestIterate[level] = (i+1) % MAX_COMBINED_ITERATES;


	// update numIterates
	if (numActiveIterates[level] < MAX_COMBINED_ITERATES) {
		numActiveIterates[level]++;
	}

	// update history array
	history[level][i] = x;
	rHistory[level][i] = r;
}

void SolverLamg::clearHistory(index level) {
	latestIterate[level] = 0;
	numActiveIterates[level] = 0;
}

void SolverLamg::minRes(index level, Vector &x, const Vector &r) {
	if (numActiveIterates[level] > 0) {
		count n = numActiveIterates[level];
		std::vector<std::pair<index, index>> AEpos;
		std::vector<double> AEvalues;

		std::vector<std::pair<index, index>> Epos;
		std::vector<double> Evalues;

		for (index i = 0; i < r.getDimension(); ++i) {
			for (index k = 0; k < n; ++k) {
				double AEvalue = r[i] - rHistory[level][k][i];
				if (std::abs(AEvalue) > 1e-9) {
					AEpos.push_back(std::make_pair(i, k));
					AEvalues.push_back(AEvalue);
				}

				double Eval = history[level][k][i] - x[i];
				if (std::abs(Eval) > 1e-9) {
					Epos.push_back(std::make_pair(i, k));
					Evalues.push_back(Eval);
				}
			}
		}


		CSRMatrix AE(r.getDimension(), n, AEpos, AEvalues);
		CSRMatrix E(r.getDimension(), n, Epos, Evalues);

		Vector alpha = smoother.relax(CSRMatrix::mTmMultiply(AE, AE), CSRMatrix::mTvMultiply(AE, r), Vector(n, 0.0), SETUP_RELAX_COARSEST_SWEEPS);
		x += E * alpha;
	}
}

} /* namespace NetworKit */
