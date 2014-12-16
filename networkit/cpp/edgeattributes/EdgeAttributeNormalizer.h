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
	EdgeAttributeNormalizer(const Graph &G, const std::vector<A> &attribute, bool invert = false, double lower = 0, double upper = 1.0) :
		G(G), input(attribute), hasOutput(false), invert(invert), lower(lower), upper(upper) {}

	void run() {
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
			offset = upper - minValue * factor;
		}

		output.resize(G.upperEdgeIdBound(), std::numeric_limits<double>::quiet_NaN());

		G.parallelForEdges([&](node u, node v, edgeid eid) {
			output[eid] = factor * input[eid] + offset;
		});

		hasOutput = true;
	}

	virtual std::vector<double> getAttribute() {
		if (!hasOutput) throw std::runtime_error("Error: Run must be called first");

		hasOutput = false;
		return std::move(output);
	}

private:
	const Graph &G;
	const std::vector<A> &input;
	std::vector<double> output;
	bool hasOutput, invert;
	double lower, upper;
};

}

#endif /* EDGEATTRIBUTENORMALIZER_H */
