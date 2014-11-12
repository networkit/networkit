#ifndef ATTRIBUTEBLENDER_H
#define ATTRIBUTEBLENDER_H

#include "../graph/Graph.h"

namespace NetworKit {

class AttributeBlender {
public:
	AttributeBlender(const Graph &G, const std::vector<double> &attribute0, const std::vector<double> &attribute1, const std::vector<bool> &selection);

	void run();

	std::vector<double> getAttribute();

private:
	const Graph &G;
	const std::vector<double> &attribute0, &attribute1;
	const std::vector<bool> &selection;
	std::vector<double> blendedAttribute;
	bool hasAttribute;
};

} // namespace NetworKit

#endif // ATTRIBUTEBLENDER_H
