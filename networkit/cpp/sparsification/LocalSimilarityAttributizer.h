/*
 * LocalSimilarityAttributizer.h
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#ifndef LOCALSIMATTRIBUTIZER_H_
#define LOCALSIMATTRIBUTIZER_H_

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

template<typename T>
struct AttributizedEdge {
	node ego;
	node alter;
	edgeid eid;
	T value;

	AttributizedEdge(node ego, node alter, edgeid eid, T v) :
			ego(ego), alter(alter), eid(eid), value(v) {
	}

	bool operator<(const AttributizedEdge<T>& other) const {
		return (value > other.value)
				|| (value == other.value && alter < other.alter);
	}

	bool operator>(const AttributizedEdge<T>& other) const {
		return (value < other.value)
				|| (value == other.value && alter > other.alter);
	}

	bool operator==(const AttributizedEdge<T>& other) const {
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
class LocalSimilarityAttributizer : public EdgeAttribute<double> {

public:

	/**
	 * Creates a new instance of the Local Sparsification algorithm.
	 */
	LocalSimilarityAttributizer(const Graph &graph, const std::vector<count>& triangles);

	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;
	const std::vector<count>& triangles;
};

}
/* namespace NetworKit */

#endif /* LOCALSIMATTRIBUTIZER_H_ */
