/*
 * EdgeAttributeNormalizer.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include "EdgeAttributeNormalizer.h"

namespace NetworKit {

template<typename A>
EdgeAttributeNormalizer<A>::EdgeAttributeNormalizer(const Graph &G, const std::vector<A> &attribute, bool invert, double lower, double upper) :
		G(G), input(attribute), hasOutput(false), invert(invert), lower(lower), upper(upper) {};

template <typename A>
void EdgeAttributeNormalizer<A>::run() {
	A minValue = std::numeric_limits< A >::max();
	A maxValue = std::numeric_limits< A >::lowest();

	G.forEdges([&](node u, node v, edgeid eid) {
		if (input[eid] < minValue) {
			minValue = input[eid];
		}

		if (input[eid] > maxValue) {
			maxValue = input[eid];
		}
	});

	double factor = (upper - lower) / (maxValue - minValue), offset = lower - minValue * factor;

	if (invert) {
		factor *= -1.0;
		offset += lower + upper;
	}

	output.resize(G.upperEdgeIdBound(), std::numeric_limits<double>::quiet_NaN());

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		output[eid] = factor * input[eid] + offset;
	});

	hasOutput = true;
}

template <typename A>
std::vector<double> EdgeAttributeNormalizer<A>::getAttribute() {
	if (!hasOutput) throw std::runtime_error("Error: Run must be called first");

	hasOutput = false;
	return std::move(output);
}

} /* namespace NetworKit */
