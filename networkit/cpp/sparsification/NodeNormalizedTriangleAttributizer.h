/*
 * NodeNormalizedTriangleAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef NODENORMALIZEDTRIANGLEATTRIBUTIZER_H
#define NODENORMALIZEDTRIANGLEATTRIBUTIZER_H

#include "../edgescores/EdgeAttribute.h"

namespace NetworKit {

class NodeNormalizedTriangleAttributizer : public EdgeAttribute<double> {

public:
	NodeNormalizedTriangleAttributizer(const Graph& graph, const std::vector<count>& triangles);
	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;
	const std::vector<count>& triangles;

};

} // namespace NetworKit

#endif // NODENORMALIZEDTRIANGLEATTRIBUTIZER_H
