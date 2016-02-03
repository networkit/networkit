/*
 * LevelElimination.cpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#include "LevelElimination.h"
#include "../../../io/LineFileReader.h"

namespace NetworKit {

count LevelElimination::lvl = 3;

LevelElimination::LevelElimination(const CSRMatrix &A, const std::vector<EliminationStage> &coarseningStages) : Level(ELIMINATION, A), coarseningStages(coarseningStages) {
	cIndexFine = std::vector<index>(A.numberOfRows());
#pragma omp parallel for
	for (index i = 0; i < cIndexFine.size(); ++i) {
		cIndexFine[i] = i;
	}

	for (index k = coarseningStages.size(); k-- > 0;) {
		for (index i = 0; i < cIndexFine.size(); ++i) {
			assert(cIndexFine[i] < coarseningStages[k].getCSet().size());
			cIndexFine[i] = coarseningStages[k].getCSet()[cIndexFine[i]];
		}
	}
}

void LevelElimination::coarseType(const Vector &xf, Vector &xc) const {
	xc = Vector(A.numberOfRows());
#pragma omp parallel for
	for (index i = 0; i < xc.getDimension(); ++i) {
		xc[i] = xf[cIndexFine[i]];
	}
}

void LevelElimination::restrict(const Vector &bf, Vector &bc, std::vector<Vector> &bStages) const {
	bStages.resize(coarseningStages.size() + 1);
	bStages[0] = bf;
	bc = bf;
	index curStage = 0;
	for (EliminationStage s : coarseningStages) {
		//Vector bOld = bStages[curStage];
		Vector bCSet;
		subVectorExtract(bCSet, bc, s.getCSet());

		Vector bFSet;
		subVectorExtract(bFSet, bc, s.getFSet());
		bc = bCSet + s.getR() * bFSet;
		bStages[curStage+1] = bc; // b = b.c + s.P^T * b.f

		curStage++;
	}
}

void LevelElimination::interpolate(const Vector &xc, Vector &xf, const std::vector<Vector> &bStages) const {
	Vector currX = xc;
	for (index k = coarseningStages.size(); k-- > 0;) {
		EliminationStage s = coarseningStages[k];
		xf = Vector(s.getN());
		Vector bFSet;
		subVectorExtract(bFSet, bStages[k], s.getFSet());

		Vector bq(bFSet.getDimension());
		const Vector &q = s.getQ();
#pragma omp parallel for
		for (index i = 0; i < bq.getDimension(); ++i) { // bq = s.q .* b.f
			bq[i] = q[i] * bFSet[i];
		}
		Vector xFSet = s.getP() * currX + bq;

		const std::vector<index> &fSet = s.getFSet();
#pragma omp parallel for
		for (index i = 0; i < xFSet.getDimension(); ++i) {
			xf[fSet[i]] = xFSet[i];
		}

		const std::vector<index> &cSet = s.getCSet();
#pragma omp parallel for
		for (index i = 0; i < currX.getDimension(); ++i) {
			xf[cSet[i]] = currX[i];
		}

		currX = xf;
	}
}

void LevelElimination::subVectorExtract(Vector &subVector, const Vector &vector, const std::vector<index> &elements) const {
	subVector = Vector(elements.size());
#pragma omp parallel for
	for (index i = 0; i < elements.size(); ++i) {
		subVector[i] = vector[elements[i]];
	}
}


} /* namespace NetworKit */
