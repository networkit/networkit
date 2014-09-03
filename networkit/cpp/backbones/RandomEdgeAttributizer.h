#ifndef RANDOMEDGEATTRIBUTIZER_H
#define RANDOMEDGEATTRIBUTIZER_H

#include "AttributeGenerator.h"

namespace NetworKit {

class RandomEdgeAttributizer : public AttributeGenerator<int, double> {
public:
	RandomEdgeAttributizer(double rneRatio = 0.8);
	
	virtual std::vector<double> getAttribute(const Graph& g, const std::vector<int> &attribute) override;
	
private:
	double rneRatio;
};

} // namespace NetworKit

#endif // RANDOMEDGEATTRIBUTIZER_H
