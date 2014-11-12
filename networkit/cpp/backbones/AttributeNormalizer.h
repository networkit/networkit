#ifndef ATTRIBUTENORMALIZER_H
#define ATTRIBUTENORMALIZER_H

#include "../graph/Graph.h"

namespace NetworKit {

template <typename A>
class AttributeNormalizer {
public:
	AttributeNormalizer(const Graph &G, const std::vector<A> &attribute, bool invert = false, double lower = 0, double upper = 1.0) :
		G(G), input(attribute), hasOutput(false), invert(invert), lower(lower), upper(upper) {};

	void run();

	std::vector<double> getAttribute() {
		if (!hasOutput) throw std::runtime_error("Error: Run must called first");

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

template <typename A>
void AttributeNormalizer<A>::run() {
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

}

#endif // ATTRIBUTENORMALIZER_H
