/*
 * MultigridHierarchyBuilder.cpp
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MultigridHierarchyBuilder.h"

namespace NetworKit {

MultigridHierarchy MultigridHierarchyBuilder::buildHierarchy(const Graph &graph) const {
	return buildHierarchy(LaplacianMatrix(graph));
}

} /* namespace NetworKit */
