/*
 * SolverLamg.cpp
 *
 *  Created on: 12.01.2015
 *      Author: Michael
 */

#include "SolverLamg.h"
#include "LAMGSettings.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "../../auxiliary/Enforce.h"
#include "../../auxiliary/Timer.h"

namespace NetworKit {

#ifndef NDEBUG
count SolverLamg::minResTime = 0;
count SolverLamg::interpolationTime = 0;
count SolverLamg::restrictionTime = 0;
count SolverLamg::coarsestSolve = 0;
#endif

SolverLamg::SolverLamg(LevelHierarchy &hierarchy, const Smoother &smoother) : hierarchy(hierarchy), smoother(smoother), bStages(hierarchy.size(), std::vector<Vector>()) {
}

void SolverLamg::solve(Vector &x, const Vector &b, LAMGSolverStatus &status) {
	bStages = std::vector<std::vector<Vector>>(hierarchy.size(), std::vector<Vector>());
	if (hierarchy.size() >= 2) {
		Vector bc = b;
		Vector xc = x;
		int finest = 0;

		if (hierarchy.getType(1) == ELIMINATION) {
#ifndef NDEBUG
			Aux::Timer t; t.start();
#endif
			hierarchy.at(1).restrict(b, bc, bStages[1]);
			if (hierarchy.at(1).getLaplacian().numberOfRows() == 1) {
				x = 0.0;
				return;
			} else {
				hierarchy.at(1).coarseType(x, xc);
				finest = 1;
			}
#ifndef NDEBUG
			t.stop();
			restrictionTime += t.elapsedMicroseconds();
#endif
		}
		solveCycle(xc, bc, finest, status);

		if (finest == 1) { // interpolate from finest == ELIMINATION level back to actual finest level
#ifndef NDEBUG
			Aux::Timer t; t.start();
#endif
			hierarchy.at(1).interpolate(xc, x, bStages[1]);
#ifndef NDEBUG
			t.stop();
			interpolationTime += t.elapsedMicroseconds();
#endif
		} else {
			x = xc;
		}
	} else {
		solveCycle(x, b, 0, status);
	}

	double residual = (b - hierarchy.at(0).getLaplacian() * x).length();
	status.residual = residual;
#ifndef NDEBUG
	DEBUG("final residual\t ", residual);
	DEBUG("minResTime: ", minResTime / 1000);
	DEBUG("interpolationTime: ", interpolationTime / 1000);
	DEBUG("restrictionTime: ", restrictionTime / 1000);
	DEBUG("coarsestSolve: ", coarsestSolve / 1000);
#endif
}

void SolverLamg::solveCycle(Vector &x, const Vector &b, int finest, LAMGSolverStatus &status) {
	Aux::Timer timer;
	timer.start();

	// data structures for iterate recombination
	history = std::vector<std::vector<Vector>>(hierarchy.size());
	rHistory = std::vector<std::vector<Vector>>(hierarchy.size());
	latestIterate = std::vector<index>(hierarchy.size(), 0);
	numActiveIterates = std::vector<count>(hierarchy.size(), 0);
	int coarsest = hierarchy.size() - 1;
	std::vector<count> numVisits(coarsest);
	std::vector<Vector> X(hierarchy.size());
	std::vector<Vector> B(hierarchy.size());

	for (index i = 0; i < hierarchy.size(); ++i) {
		history[i] = std::vector<Vector>(MAX_COMBINED_ITERATES, Vector(hierarchy.at(i).getNumberOfNodes()));
		rHistory[i] = std::vector<Vector>(MAX_COMBINED_ITERATES, Vector(hierarchy.at(i).getNumberOfNodes()));
	}

	Vector r = b - hierarchy.at(finest).getLaplacian() * x;
	double residual = r.length();
	double finalResidual = residual * status.desiredResidualReduction;
	double bestResidual = std::numeric_limits<double>::max();

	count iterations = 0;
	status.residualHistory.emplace_back(residual);
	count noResReduction = 0;
	while (residual > finalResidual && noResReduction < 5 && iterations < status.maxIters && timer.elapsedMilliseconds() <= status.maxConvergenceTime ) {
#ifndef NDEBUG
		DEBUG("iter ", iterations, " r=", residual);
#endif
		cycle(x, b, finest, coarsest, numVisits, X, B, status);
		r = b - hierarchy.at(finest).getLaplacian() * x;
		residual = r.length();
		status.residualHistory.emplace_back(residual);
		if (residual < bestResidual) {
			noResReduction = 0;
			bestResidual = residual;
		} else {
			++noResReduction;
		}
		iterations++;
	}

	timer.stop();

	status.numIters = iterations;
	status.residual = r.length();
	status.converged = r.length() <= finalResidual;
#ifndef NDEBUG
	DEBUG("nIter\t ", iterations);
#endif
}

void SolverLamg::cycle(Vector &x, const Vector &b, int finest, int coarsest, std::vector<count> &numVisits, std::vector<Vector> &X, std::vector<Vector> &B, const LAMGSolverStatus &status) {
	std::fill(numVisits.begin(), numVisits.end(), 0);
	X[finest] = x;
	B[finest] = b;

#ifndef NDEBUG
	Aux::Timer t;
#endif

	int currLvl = finest;
	int nextLvl = finest;
	double maxVisits = 0.0;

	saveIterate(currLvl, X[currLvl], B[currLvl] - hierarchy.at(currLvl).getLaplacian() * X[currLvl]);
	while (true) {
		if (currLvl == coarsest) {
#ifndef NDEBUG
			t.start();
#endif
			nextLvl = currLvl - 1;
			if (currLvl == finest) { // finest level
				X[currLvl] = smoother.relax(hierarchy.at(currLvl).getLaplacian(), B[currLvl], X[currLvl], status.numPreSmoothIters);
			} else {
				Vector bCoarse(B[currLvl].getDimension()+1, 0.0);
				for (index i = 0; i < B[currLvl].getDimension(); ++i) {
					bCoarse[i] = B[currLvl][i];
				}

				Vector xCoarse = DenseMatrix::LUSolve(hierarchy.getCoarseMatrix(), bCoarse);
				for (index i = 0; i < X[currLvl].getDimension(); ++i) {
					X[currLvl][i] = xCoarse[i];
				}
			}
#ifndef NDEBUG
			t.stop();
			coarsestSolve += t.elapsedMicroseconds();
#endif
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
#ifndef NDEBUG
			t.start();
#endif
			numVisits[currLvl]++;

			if (hierarchy.getType(nextLvl) != ELIMINATION) {
				X[currLvl] = smoother.relax(hierarchy.at(currLvl).getLaplacian(), B[currLvl], X[currLvl], status.numPreSmoothIters);
			}

			if (hierarchy.getType(nextLvl) == ELIMINATION) {
				hierarchy.at(nextLvl).restrict(B[currLvl], B[nextLvl], bStages[nextLvl]);
			} else {
				hierarchy.at(nextLvl).restrict(B[currLvl] - hierarchy.at(currLvl).getLaplacian() * X[currLvl], B[nextLvl]);
			}

			hierarchy.at(nextLvl).coarseType(X[currLvl], X[nextLvl]);

			clearHistory(nextLvl);
#ifndef NDEBUG
			t.stop();
			restrictionTime += t.elapsedMicroseconds();
#endif
		} else { // postProcess
			if (currLvl == coarsest || hierarchy.getType(currLvl+1) != ELIMINATION) {
				minRes(currLvl, X[currLvl], B[currLvl] - hierarchy.at(currLvl).getLaplacian() * X[currLvl]);
			}

#ifndef NDEBUG
			t.start();
#endif

			if (nextLvl > finest) {
				saveIterate(nextLvl, X[nextLvl], B[nextLvl] - hierarchy.at(nextLvl).getLaplacian() * X[nextLvl]);
			}


			if (hierarchy.getType(currLvl) == ELIMINATION) {
				hierarchy.at(currLvl).interpolate(X[currLvl], X[nextLvl], bStages[currLvl]);
			} else {
				Vector xf = X[nextLvl];
				hierarchy.at(currLvl).interpolate(X[currLvl], xf);
				X[nextLvl] += xf;
			}

			if (hierarchy.getType(currLvl) != ELIMINATION) {
				X[nextLvl] = smoother.relax(hierarchy.at(nextLvl).getLaplacian(), B[nextLvl], X[nextLvl], status.numPostSmoothIters);
			}

#ifndef NDEBUG
			t.stop();
			interpolationTime += t.elapsedMicroseconds();
#endif
		}

		currLvl = nextLvl;
	} // while

	// post-cycle finest
	if ((int64_t) hierarchy.size() > finest + 1 && hierarchy.getType(finest+1) != ELIMINATION) { // do an iterate recombination on calculated solutions
		minRes(finest, X[finest], B[finest] - hierarchy.at(finest).getLaplacian() * X[finest]);
	}


	X[finest] -= X[finest].mean();
	x = X[finest];
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

		std::vector<index> ARowIdx(r.getDimension()+1);
		std::vector<index> ERowIdx(r.getDimension()+1);

#pragma omp parallel for
		for (index i = 0; i < r.getDimension(); ++i) {
			for (index k = 0; k < n; ++k) {
				double AEvalue = r[i] - rHistory[level][k][i];
				if (std::abs(AEvalue) > 1e-9) {
					++ARowIdx[i+1];
				}

				double Eval = history[level][k][i] - x[i];
				if (std::abs(Eval) > 1e-9) {
					++ERowIdx[i+1];
				}
			}
		}

		for (index i = 0; i < r.getDimension(); ++i) {
			ARowIdx[i+1] += ARowIdx[i];
			ERowIdx[i+1] += ERowIdx[i];
		}

		std::vector<index> AColumnIdx(ARowIdx[r.getDimension()]);
		std::vector<double> ANonZeros(ARowIdx[r.getDimension()]);

		std::vector<index> EColumnIdx(ERowIdx[r.getDimension()]);
		std::vector<double> ENonZeros(ERowIdx[r.getDimension()]);

#pragma omp parallel for
		for (index i = 0; i < r.getDimension(); ++i) {
			for (index k = 0, aIdx = ARowIdx[i], eIdx = ERowIdx[i]; k < n; ++k) {
				double AEvalue = r[i] - rHistory[level][k][i];
				if (std::abs(AEvalue) > 1e-9) {
					AColumnIdx[aIdx] = k;
					ANonZeros[aIdx] = AEvalue;
					++aIdx;
				}

				double Eval = history[level][k][i] - x[i];
				if (std::abs(Eval) > 1e-9) {
					EColumnIdx[eIdx] = k;
					ENonZeros[eIdx] = Eval;
					++eIdx;
				}
			}
		}

		CSRMatrix AE(r.getDimension(), n, ARowIdx, AColumnIdx, ANonZeros, true);
		CSRMatrix E(r.getDimension(), n, ERowIdx, EColumnIdx, ENonZeros, true);
#ifndef NDEBUG
	Aux::Timer t;
	t.start();
#endif

		Vector alpha = smoother.relax(CSRMatrix::mTmMultiply(AE, AE), CSRMatrix::mTvMultiply(AE, r), Vector(n, 0.0), 10);
		x += E * alpha;

#ifndef NDEBUG
	t.stop();
	minResTime += t.elapsedMicroseconds();
#endif
	}

}

} /* namespace NetworKit */
