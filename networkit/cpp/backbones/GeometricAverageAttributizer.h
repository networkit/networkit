#ifndef GEOMETRICAVERAGEATTRIBUTIZER_H
#define GEOMETRICAVERAGEATTRIBUTIZER_H

#include "AttributeGenerator.h"

namespace NetworKit {

class GeometricAverageAttributizer : public AttributeGenerator<double, double> {
public:
	virtual std::vector<double> getAttribute(const Graph& g, const std::vector<double>& attribute) override;
};

} // namespace NetworKit

#endif // GEOMETRICAVERAGEATTRIBUTIZER_H
