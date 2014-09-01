#ifndef CHANCECORRECTEDTRIANGLEATTRIBUTIZER_H
#define CHANCECORRECTEDTRIANGLEATTRIBUTIZER_H

#include "AttributeGenerator.h"

namespace NetworKit {

class ChanceCorrectedTriangleAttributizer : public AttributeGenerator<int, double> {
public:
	virtual std::vector<double> getAttribute(const Graph& graph, const std::vector<int>& triangles) override;

};

}
#endif // CHANCECORRECTEDTRIANGLEATTRIBUTIZER_H
