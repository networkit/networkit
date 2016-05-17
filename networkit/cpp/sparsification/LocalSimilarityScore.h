/*
 * LocalSimilarityScore.h
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#ifndef LOCALSIMSCORE_H_
#define LOCALSIMSCORE_H_

#include "../edgescores/EdgeScore.h"

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
class LocalSimilarityScore : public EdgeScore<double> {

public:

	LocalSimilarityScore(const Graph& G, const std::vector<count>& triangles);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

private:
	const std::vector<count>& triangles;
};

}
/* namespace NetworKit */

#endif /* LOCALSIMSCORE_H_ */
