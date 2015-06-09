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
	std::vector<Matrix> laplacianMatrices;
	std::vector<Matrix> interpolationMatrices;

public:
	MultigridHierarchy() {}
	MultigridHierarchy(const Matrix &fineMatrix);
	virtual ~MultigridHierarchy() {}

	void addLevel(const Matrix &coarseMatrix, const Matrix &interpolationMatrix);

	const Matrix& getLaplacian(const index level);
	const Matrix& getInterpolationMatrix(const index level);

	virtual count getNumPreSmooth(const index level) const;
	virtual count getNumPostSmooth(const index level) const;
	virtual float getNumMultigridCycles(const index level) const;

	virtual Vector createInitialCoarseResult(const index level, const Vector &xFine) const;
	virtual Vector restriction(const index level, const Vector &xFine, const Vector &bFine) const;
	virtual Vector prolongation(const index level, const Vector &xCoarse, const Vector &xFine, const Vector &bFine) const;

	count getNumLevels() const;

};

} /* namespace NetworKit */

#endif /* MULTIGRIDHIERARCHY_H_ */
