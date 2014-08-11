#ifndef ATTRIBUTEASWEIGHT_H
#define ATTRIBUTEASWEIGHT_H

#include "../graph/Graph.h"

namespace NetworKit {

class AttributeAsWeight {
public:
	AttributeAsWeight(bool squared = false, edgeweight offset = 1, edgeweight factor = 1);
	virtual ~AttributeAsWeight() = default;
	virtual Graph calculate(Graph &g, const std::vector<double> &attribute);
private:
	bool squared;
	edgeweight offset;
	edgeweight factor;
};

} // namespace NetworKit

#endif // ATTRIBUTEASWEIGHT_H
