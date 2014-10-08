#ifndef CHIBANISHIZEKIQUADRANGLECOUNTER_H
#define CHIBANISHIZEKIQUADRANGLECOUNTER_H

#include "AttributeGenerator.h"

namespace NetworKit {

class ChibaNishizekiQuadrangleCounter : public AttributeGenerator<int, int> {
public:
	virtual std::vector<int> getAttribute(const Graph& graph, const std::vector<int>& attribute) override;
};

} // namespace NetworKit

#endif // CHIBANISHIZEKIQUADRANGLECOUNTER_H
