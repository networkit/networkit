#ifndef LINEARIZEATTRIBUTE_H
#define LINEARIZEATTRIBUTE_H

#include "AttributeGenerator.h"

namespace NetworKit {

class LinearizeAttribute : public AttributeGenerator<double, double> {
public:
	virtual std::vector<double> getAttribute(const NetworKit::Graph &g, const std::vector< double > &attribute) override;

};

} // namespace NetworKit

#endif // LINEARIZEATTRIBUTE_H
