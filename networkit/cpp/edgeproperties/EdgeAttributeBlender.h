/*
 * EdgeAttributeBlender.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGEATTRIBUTEBLENDER_H
#define EDGEATTRIBUTEBLENDER_H

#include "../graph/Graph.h"

namespace NetworKit {

class EdgeAttributeBlender {

public:
	EdgeAttributeBlender(const Graph &G, const std::vector<double> &attribute0, const std::vector<double> &attribute1, const std::vector<bool> &selection);

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

#endif // EDGEATTRIBUTEBLENDER_H
