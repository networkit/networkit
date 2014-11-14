/*
 * MultigridHierarchyBuilder.h
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MULTIGRIDHIERARCHYBUILDER_H_
#define MULTIGRIDHIERARCHYBUILDER_H_

#include "MultigridHierarchy.h"
#include "../algebraic/LaplacianMatrix.h"

namespace NetworKit {

class MultigridHierarchyBuilder {
public:
	MultigridHierarchyBuilder() {}
	virtual ~MultigridHierarchyBuilder() {}

	virtual MultigridHierarchy buildHierarchy(const Matrix &matrix) const = 0;
	virtual MultigridHierarchy buildHierarchy(const Graph &graph) const;
};

} /* namespace NetworKit */

#endif /* MULTIGRIDHIERARCHYBUILDER_H_ */
