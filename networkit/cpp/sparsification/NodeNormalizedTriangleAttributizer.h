/*
 * NodeNormalizedTriangleAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef NODENORMALIZEDTRIANGLEATTRIBUTIZER_H
#define NODENORMALIZEDTRIANGLEATTRIBUTIZER_H

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

class NodeNormalizedTriangleAttributizer : public EdgeAttribute<double> {

public:
	NodeNormalizedTriangleAttributizer(const Graph& graph, const std::vector<int>& triangles);
	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;
	const std::vector<int>& triangles;

};

} // namespace NetworKit

#endif // NODENORMALIZEDTRIANGLEATTRIBUTIZER_H
