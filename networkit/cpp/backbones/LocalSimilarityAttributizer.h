/*
 * LocalSimilarityAttributizer.h
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#ifndef LOCALSIMATTRIBUTIZER_H_
#define LOCALSIMATTRIBUTIZER_H_

#include "BackboneCalculator.h"
#include "gtest/gtest_prod.h"

namespace NetworKit {

struct AttributizedEdge {
	node ego;
	node alter;
	edgeid eid;
	double value;

	AttributizedEdge(node ego, node alter, edgeid eid, double v) :
			ego(ego), alter(alter), eid(eid), value(v) {
	}

	bool operator<(const AttributizedEdge& other) const {
		return (value > other.value)
				|| (value == other.value && alter < other.alter);
	}

	bool operator>(const AttributizedEdge& other) const {
		return (value < other.value)
				|| (value == other.value && alter > other.alter);
	}

	bool operator==(const AttributizedEdge& other) const {
		return ego == other.ego && alter == other.alter
				&& value == other.value;
	}
};

struct greater {
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

/**
 * Implementation of the Local Sparsification Algorithm by Sataluri et al.
 */
class LocalSimilarityAttributizer : public AttributeGenerator<int, double> {

public:

	/**
	 * Creates a new instance of the Local Sparsification algorithm.
	 */
	LocalSimilarityAttributizer();

	~LocalSimilarityAttributizer() = default;

	std::vector<double> getAttribute(const Graph& graph, const std::vector<int>& attribute);

private:
	//Private helper functions
	double getSimilarity(const Graph& graph, node i, node j);

	FRIEND_TEST(LocalSimilarityGTest, testSimilarityCalculation);
};

}
/* namespace NetworKit */
#endif /* LOCALSIMATTRIBUTIZER_H_ */
