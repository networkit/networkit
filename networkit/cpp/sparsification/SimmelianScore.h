/*
 * SimmelianAttributizer.h
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#ifndef SIMMELIANATTRIBUTIZER_H_
#define SIMMELIANATTRIBUTIZER_H_

#include "../edgescores/EdgeScore.h"
#include <set>

namespace NetworKit {

/**
 * A directed edge with a simmelianness int value.
 * An ordering is defined on these ties: We order first by
 * simmelianness (descending) and then by the id of alter (ascending).
 */
struct RankedEdge {
	node ego;
	node alter;
	count simmelianness; 	//The number of triangles the edge is embedded in.
	count rank; 			//The rank within the ranked neighborhood.

	RankedEdge(node ego, node alter, count s) :
			ego(ego), alter(alter), simmelianness(s), rank(0) {
	}

	RankedEdge(node ego, node alter, count s, count r) :
			ego(ego), alter(alter), simmelianness(s), rank(r) {
	}

	bool operator<(const RankedEdge& other) const {
		return (simmelianness > other.simmelianness)
				|| (simmelianness == other.simmelianness && alter > other.alter);
	}

	bool operator>(const RankedEdge& other) const {
		return (simmelianness < other.simmelianness)
				|| (simmelianness == other.simmelianness && alter < other.alter);
	}

	bool operator==(const RankedEdge& other) const {
		return ego == other.ego && simmelianness == other.simmelianness
				&& alter == other.alter && rank == other.rank;
	}
};

typedef std::vector<RankedEdge> RankedNeighbors;

/**
 * Represents the result of the comparison of two ranked neighborhood lists, namely
 * the overlap of the top-k neighbors and the maximum jaccard index.
 */
struct Redundancy {
	count overlap;
	double jaccard;

	Redundancy(count o, double j) : overlap(o), jaccard(j) { }
};

/**
 * Abstract base class for the two variants of Simmelian backbones (OverlapFilter, JaccardFilter).
 */
class SimmelianScore : public EdgeScore<double> {

public:

	SimmelianScore(const Graph& graph, const std::vector<count>& attribute);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;

	std::vector<RankedNeighbors> getRankedNeighborhood(const Graph& g, const std::vector<count>& triangles);

	Redundancy getOverlap(
			const node& ego,
			const node& alter,
			const std::vector<RankedNeighbors>& neighbors,
			const count& maxRank);

	void matchNeighbors(
			const node& ego,
			const node& alter,
			const bool& reciprocityCheck,
			std::vector<RankedEdge>::const_iterator& egoIt,
			const RankedNeighbors& egoNeighbors,
			std::set<node>& egoNeighborsUnmatched,
			std::set<node>& alterNeighborsUnmatched,
			const count& rank,
			count& overlap);

protected:
	const std::vector<count>& triangles;

};

} /* namespace NetworKit */

#endif /* SIMMELIANSCORE_H_ */
