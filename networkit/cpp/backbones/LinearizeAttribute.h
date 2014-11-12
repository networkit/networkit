#ifndef LINEARIZEATTRIBUTE_H
#define LINEARIZEATTRIBUTE_H

#include "AttributeGenerator.h"

namespace NetworKit {

class LinearizeAttribute : public AttributeGenerator<double, double> {
public:
	LinearizeAttribute(bool inverse = false);
	virtual std::vector<double> getAttribute(const Graph& g, const std::vector<double>& attribute) override;
private:
	bool inverse;
};

} // namespace NetworKit

#endif // LINEARIZEATTRIBUTE_H
