/*
 * LAMGHierarchy.h
 *
 *  Created on: 14.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LAMGHIERARCHY_H_
#define LAMGHIERARCHY_H_

#include "MultigridHierarchy.h"
#include <unordered_map>

namespace NetworKit {

class LAMGHierarchy : public MultigridHierarchy {
public:
	enum Type{ELIMINATION, AGGREGATION};

	LAMGHierarchy() {}
	LAMGHierarchy(const Matrix &fineMatrix);

	void addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, Type type);
	void addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, const Matrix &qMatrix);
	using MultigridHierarchy::addLevel;

	count getNumPreSmooth(const index level) const override;
	count getNumPostSmooth(const index level) const override;
	float getNumMultigridCycles(const index level) const override;

	Vector createInitialCoarseResult(const index level, const Vector &xFine) const;
	Vector restriction(const index level, const Vector &xFine, const Vector &bFine) const override;
	Vector prolongation(const index level, const Vector &xCoarse, const Vector &xFine, const Vector &bFine) const override;

	Type getType(const index level) const;
	count nnzOfMatrix(const index level) const;

private:
	std::vector<Type> type;
	std::unordered_map<index, Matrix> qMatrix;
};

} /* namespace NetworKit */

#endif /* LAMGHIERARCHY_H_ */
