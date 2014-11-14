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
private:
	std::vector<Matrix> laplacians;
	std::vector<Matrix> interpolationMatrices;

public:
	MultigridHierarchy() {}

	void addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix);

	const Matrix& getLaplacian(const index level);
	const Matrix& getInterpolationMatrix(const index level);

	count getNumLevels() const;

};

} /* namespace NetworKit */

#endif /* MULTIGRIDHIERARCHY_H_ */
