/*
 * EdgeAttributeBlender.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGEATTRIBUTEBLENDER_H
#define EDGEATTRIBUTEBLENDER_H

#include "../graph/Graph.h"
#include "EdgeAttribute.h"

namespace NetworKit {

class EdgeAttributeBlender : public EdgeAttribute<double> {

public:

	EdgeAttributeBlender(const Graph &G, const std::vector<double> &attribute0, const std::vector<double> &attribute1, const std::vector<bool> &selection);

	void run();

	virtual std::vector<double> getAttribute() override;

private:
	const Graph &G;
	const std::vector<double> &attribute0, &attribute1;
	const std::vector<bool> &selection;
	std::vector<double> blendedAttribute;
	bool hasAttribute;
};

} // namespace NetworKit

#endif // EDGEATTRIBUTEBLENDER_H
