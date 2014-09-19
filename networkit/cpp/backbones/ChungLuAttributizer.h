#ifndef CHUNGLUATTRIBUTIZER_H
#define CHUNGLUATTRIBUTIZER_H

#include "AttributeGenerator.h"

namespace NetworKit {

class ChungLuAttributizer : public AttributeGenerator<int, double> {
public:

	virtual std::vector<double> getAttribute(const Graph& g, const std::vector<int>& attribute) override;

};

} // namespace NetworKit

#endif // CHUNGLUATTRIBUTIZER_H
