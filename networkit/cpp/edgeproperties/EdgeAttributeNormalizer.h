/*
 * EdgeAttributeNormalizer.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGEATTRIBUTENORMALIZER_H
#define EDGEATTRIBUTENORMALIZER_H

#include "../graph/Graph.h"
#include "EdgeAttribute.h"

namespace NetworKit {

template <typename A>
class EdgeAttributeNormalizer : public EdgeAttribute<double> {

public:
	EdgeAttributeNormalizer(const Graph &G, const std::vector<A> &attribute, bool invert = false, double lower = 0, double upper = 1.0);
	void run();
	virtual std::vector<double> getAttribute() override;

private:
	const Graph &G;
	const std::vector<A> &input;
	std::vector<double> output;
	bool hasOutput, invert;
	double lower, upper;
};

}

#endif /* EDGEATTRIBUTENORMALIZER_H */
