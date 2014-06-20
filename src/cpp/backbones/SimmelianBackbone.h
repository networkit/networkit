/*
 * SimmelianBackbone.h
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#ifndef SIMMELIANBACKBONE_H_
#define SIMMELIANBACKBONE_H_

#include "BackboneCalculator.h"
#include "TriangleCounter.h"
#include "gtest/gtest_prod.h"
#include <set>

namespace NetworKit {

/**
 * A directed edge with a simmelianness int value.
 * An ordering is defined on these ties by considering
 * first the simmelianness and then the index of alter.
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
				|| (simmelianness == other.simmelianness && alter < other.alter);
	}

	bool operator>(const RankedEdge& other) const {
		return (simmelianness < other.simmelianness)
				|| (simmelianness == other.simmelianness && alter > other.alter);
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
 * Calculates the simmelian backbone for a given input graph.
 */
class SimmelianBackbone : public BackboneCalculator {

public:

	/**
	 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
	 * @param g 			the graph to calculate the backbone for
	 * @param maxRank 		the maximum rank that is considered for overlap calculation
	 * @param minOverlap	a minimum number of common top-maxRank neighbors
	 */
	SimmelianBackbone(const Graph& g, count maxRank, count minOverlap);

	/**
	 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
	 * @param g 			the graph to calculate the backbone for
	 * @param jaccardTresh	the minimum best prefix jaccard coefficient for an edge to be contained in the backbone
	 */
	SimmelianBackbone(const Graph& g, double jaccardThreshold);

	Graph calculate();

private:
	//Calculation parameters
	const Graph& graph;
	bool parameterized;
	double jaccardTreshold;
	count maxRank;
	count minOverlap;

	std::vector<RankedNeighbors> getRankedNeighborhood(const Graph& g, edgeCountMap& triangles);

	Redundancy getOverlap(const RankedNeighbors& first, const RankedNeighbors& second, const count& maxRank);

	void matchNeighbors(std::vector<RankedEdge>::const_iterator& egoIt,
			const RankedNeighbors& egoNeighbors,
			std::set<node>& egoNeighborsUnmatched,
			std::set<node>& alterNeighborsUnmatched,
			const count& rank,
			count& overlap);

	FRIEND_TEST(SimmelianBackboneGTest, testOverlapCounting);
	FRIEND_TEST(SimmelianBackboneGTest, testRankedNeighborhood);
	FRIEND_TEST(SimmelianBackboneGTest, testRankedNeighborhoodSkippedRanks);
	FRIEND_TEST(SimmelianBackboneGTest, testOverlapFiltering);
};

}
/* namespace NetworKit */
#endif /* SIMMELIANBACKBONE_H_ */
