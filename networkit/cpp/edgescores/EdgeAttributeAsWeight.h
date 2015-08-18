/*
 * EdgeAttributeAsWeight.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef ATTRIBUTEASWEIGHT_H
#define ATTRIBUTEASWEIGHT_H

#include "../graph/Graph.h"

namespace NetworKit {

class EdgeAttributeAsWeight {

public:
	EdgeAttributeAsWeight(const Graph& graph, const std::vector<double>& attribute, bool squared = false, edgeweight offset = 1, edgeweight factor = 1);
	virtual ~EdgeAttributeAsWeight() = default;
	virtual Graph calculate();

private:
	const Graph& graph;
	const std::vector<double>& attribute;
	bool squared;
	edgeweight offset;
	edgeweight factor;

};

} // namespace NetworKit

#endif // ATTRIBUTEASWEIGHT_H
