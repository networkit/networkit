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

namespace NetworKit {

/**
 * A directed edge with a simmelianness int value.
 * An ordering is defined on these ties by considering
 * first the simmelianness and then the index of alter.
 */
struct RankedEdge {
	node ego;
	node alter;
	int simmelianness; 	//The number of triangles the edge is embedded in.
	int rank; 			//The rank within the ranked neighborhood.

	RankedEdge (node ego, node alter, int s) : ego(ego), alter(alter), simmelianness(s), rank(0) {
	}

	RankedEdge (node ego, node alter, int s, int r) : ego(ego), alter(alter), simmelianness(s), rank(r) {
	}

	bool operator<(const RankedEdge& other) const {
		return (simmelianness > other.simmelianness) || (simmelianness == other.simmelianness && alter < other.alter);
	}

	bool operator>(const RankedEdge& other) const {
		return (simmelianness < other.simmelianness) || (simmelianness == other.simmelianness && alter > other.alter);
	}

	bool operator==(const RankedEdge& other) const {
		return simmelianness == other.simmelianness && alter == other.alter; //TODO: compare ego as well?
	}
};

typedef std::vector<RankedEdge> RankedNeighbors;

/** 
 * Calculates the simmelian backbone for a given input graph.
 */
class SimmelianBackbone : public BackboneCalculator {

public:
	/**
	 * Calculates the simmelian backbone for the given graph.
	 * @param g 			the graph to calculate the backbone for
	 * @param maxRank 		the maximum rank that is considered for overlap calculation
	 * @param minOverlap	a minimum number of common top-maxRank neighbors
	 */
	Graph calculate(const Graph& g, const int& maxRank, const int& minOverlap);

	count test();

private:
	FRIEND_TEST(SimmelianBackboneGTest, testOverlapCounting);
	FRIEND_TEST(SimmelianBackboneGTest, testRankedNeighborhood);
	FRIEND_TEST(SimmelianBackboneGTest, testOverlapFiltering);
	std::vector<RankedNeighbors> getRankedNeighborhood(const Graph& g, edgeCountMap& triangles);
	int getOverlap(const RankedNeighbors& first, const RankedNeighbors& second, const int& maxRank);
	void matchNeighbors(std::vector<RankedEdge>::const_iterator& egoIt,
						const RankedNeighbors& egoNeighbors,
						std::set<node>& egoNeighborsUnmatched,
						std::set<node>& alterNeighborsUnmatched,
						const int& rank,
						int& overlap);
};

} /* namespace NetworKit */
#endif /* SIMMELIANBACKBONE_H_ */
