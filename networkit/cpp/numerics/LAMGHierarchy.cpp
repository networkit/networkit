/*
 * LAMGHierarchy.cpp
 *
 *  Created on: 14.11.2014
 *      Author: Michael
 */

#include "LAMGHierarchy.h"

namespace NetworKit {

LAMGHierarchy::LAMGHierarchy(const Matrix &fineMatrix) : MultigridHierarchy(fineMatrix) {
}

void LAMGHierarchy::addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, Type type) {
	addLevel(laplacian, interpolationMatrix);
	this->type.push_back(type);
}

void LAMGHierarchy::addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, const Matrix &qMatrix) {
	this->qMatrix.insert(std::make_pair(getNumLevels() - 1, qMatrix));
	addLevel(laplacian, interpolationMatrix);
	type.push_back(ELIMINATION);
}

LAMGHierarchy::Type LAMGHierarchy::getType(const index level) const {
	assert(level < getNumLevels());
	return type[level];
}

count LAMGHierarchy::nnzOfMatrix(const index level) const {
	assert(level < getNumLevels());
	return laplacianMatrices[level].nnz();
}

float LAMGHierarchy::getNumMultigridCycles(const index level) const {
	assert(level < getNumLevels() - 1);
	float cycles = 1.0;
	if (level > 0 && type[level + 1] != ELIMINATION) {
		count nnzFinest = nnzOfMatrix(0);
		count nnzCurrent = nnzOfMatrix(level);
		if (nnzCurrent > 0.1 * nnzFinest) {
			cycles = 1.5;
		} else {
			count nnzNext = nnzOfMatrix(level + 1);
			cycles = std::min(2.0, 0.7 * (double) nnzNext / (double) nnzCurrent);
		}
	}

	return cycles;
}

count LAMGHierarchy::getNumPreSmooth(const index level) const {
	assert(level < getNumLevels() - 1);
	if (type[level] == ELIMINATION) {
		return 0;
	} else {
		return 1;
	}
}

count LAMGHierarchy::getNumPostSmooth(const index level) const {
	assert(level < getNumLevels() - 1);
	if (type[level] == ELIMINATION) {
		return 0;
	} else {
		return 2;
	}
}

Vector LAMGHierarchy::restriction(const index level, const Vector &xFine, const Vector &bFine) const {
	assert(level < getNumLevels());
	Vector bCoarse;

	switch(type[level]) {
		case ELIMINATION:
			bCoarse = Matrix::mTvMultiply(interpolationMatrices[level], bFine); // coarse rhs
			break;
		case AGGREGATION:
			bCoarse = Matrix::mTvMultiply(interpolationMatrices[level], bFine - laplacianMatrices[level] * xFine); // coarse residual
			break;
		default:
			break;
	}

	return bCoarse;
}

Vector LAMGHierarchy::createInitialCoarseResult(const index level, const Vector &xFine) const {
	assert(level < getNumLevels());
	if (type[level] == ELIMINATION) {
		Vector coarseResult(interpolationMatrices[level].numberOfColumns());
		count j = 0;
		for (index i = 0; i < xFine.getDimension(); ++i) {
			if (interpolationMatrices[level].nnzInRow(i) == 1) {
				coarseResult[j] = xFine[i];
				j++;
			}
		}

//		INFO("Dim=", coarseResult.getDimension(), ", j=", j);

		return coarseResult;
	} else {
		return Vector(interpolationMatrices[level].numberOfColumns(), 0.0);
	}
}

Vector LAMGHierarchy::prolongation(const index level, const Vector &xCoarse, const Vector &xFine, const Vector &bFine) const {
	assert(level < getNumLevels());

	if (type[level] == ELIMINATION) {
//		Vector bQ(bFine.getDimension());
//		Vector q = qVector.at(level);
//		for (index i = 0; i < bFine.getDimension(); ++i) {
//			bQ[i] = q[i] * bFine[i];
//		}


		return interpolationMatrices[level] * xCoarse + qMatrix.at(level) * bFine;
	} else {
		return xFine + interpolationMatrices[level] * xCoarse;
	}
}

} /* namespace NetworKit */
