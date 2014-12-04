/*
 * LAMGHierarchy.h
 *
 *  Created on: 14.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LAMGHIERARCHY_H_
#define LAMGHIERARCHY_H_

#include "MultigridHierarchy.h"

namespace NetworKit {

class LAMGHierarchy : public MultigridHierarchy {
public:
	enum Type{ELIMINATION, AGGREGATION};

	LAMGHierarchy() {}
	LAMGHierarchy(const Matrix &fineMatrix);

	void addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, Type type);
	using MultigridHierarchy::addLevel;

	count getNumPreSmooth(const index level) const;
	count getNumPostSmooth(const index level) const;
	count getNumMultigridCycles(const index level) const;

	Vector restriction(const index level, const Vector &currentApproximation, const Vector &rhs) const;
	Vector prolongation(const index level, const Vector &coarseApproximation, const Vector &fineApproximation) const;

	Type getType(const index level) const;
	count nnzOfMatrix(const index level) const;

private:
	std::vector<Type> type;
};

} /* namespace NetworKit */

#endif /* LAMGHIERARCHY_H_ */
