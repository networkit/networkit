#ifndef NODENORMALIZEDTRIANGLEATTRIBUTIZER_H
#define NODENORMALIZEDTRIANGLEATTRIBUTIZER_H

#include "AttributeGenerator.h"

namespace NetworKit {

class NodeNormalizedTriangleAttributizer : public AttributeGenerator<int, double> {
public:
	virtual std::vector<double> getAttribute(const Graph& graph, const std::vector<int>& triangles) override;
};

} // namespace NetworKit

#endif // NODENORMALIZEDTRIANGLEATTRIBUTIZER_H
