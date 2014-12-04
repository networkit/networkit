/*
 * MultigridHierarchy.h
 *
 *  Created on: 30.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MULTIGRIDHIERARCHY_H_
#define MULTIGRIDHIERARCHY_H_

#include "../algebraic/Matrix.h"

namespace NetworKit {

class MultigridHierarchy {
protected:
	std::vector<Matrix> coarseMatrices;
	std::vector<Matrix> interpolationMatrices;

public:
	MultigridHierarchy() {}
	MultigridHierarchy(const Matrix &fineMatrix);

	void addLevel(const Matrix &coarseMatrix, const Matrix &interpolationMatrix);

	const Matrix& getLaplacian(const index level);
	const Matrix& getInterpolationMatrix(const index level);

	virtual count getNumPreSmooth(const index level) const;
	virtual count getNumPostSmooth(const index level) const;
	virtual count getNumMultigridCycles(const index level) const;

	Vector restriction(const index level, const Vector &currentApproximation, const Vector &rhs) const;
	Vector prolongation(const index level, const Vector &coarseApproximation, const Vector &fineApproximation) const;

	count getNumLevels() const;

};

} /* namespace NetworKit */

#endif /* MULTIGRIDHIERARCHY_H_ */
