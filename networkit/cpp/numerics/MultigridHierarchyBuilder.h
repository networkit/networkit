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

template<class H>
class MultigridHierarchyBuilder {
public:
	MultigridHierarchyBuilder() {}
	virtual ~MultigridHierarchyBuilder() {}

	virtual H buildHierarchy(const Matrix &matrix) const = 0;
	virtual H buildHierarchy(const Graph &graph) const;
};

template <typename H>
H MultigridHierarchyBuilder<H>::buildHierarchy(const Graph &graph) const {
	return buildHierarchy(LaplacianMatrix(graph));
}

} /* namespace NetworKit */

#endif /* MULTIGRIDHIERARCHYBUILDER_H_ */
