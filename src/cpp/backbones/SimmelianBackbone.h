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

namespace NetworKit {

/**
 * An undirected edge with a simmelianness int value.
 * An ordering is defined on these ties by considering
 * first the simmelianness and then the edge.
 */
struct SimmelianTie {
	uEdge edge;
	int simmelianness; 	//The number of triangles the edge is embedded in.
	int rank; 			//The rank within the ranked neighborhood.

	SimmelianTie (uEdge e, int s) : edge(e), simmelianness(s) {
	}

	bool operator<(const SimmelianTie& other) const {
		return (simmelianness < other.simmelianness) || (simmelianness == other.simmelianness && edge < other.edge);
	}

	bool operator>(const SimmelianTie& other) const {
		return (simmelianness > other.simmelianness) || (simmelianness == other.simmelianness && edge > other.edge);
	}

	bool operator==(const SimmelianTie& other) const {
		return simmelianness == other.simmelianness && edge == other.edge;
	}
};

typedef std::vector<SimmelianTie> RankedNeighbors;

/** 
 * Calculates the simmelian backbone for a given input graph.
 */
class SimmelianBackbone : public BackboneCalculator {

public:
	/**
	 * Calculates the simmelian backbone for the given graph.
	 */
	virtual Graph calculate(const Graph& g);

	virtual count test();

private:
	std::vector<RankedNeighbors> getRankedNeighborhood(const Graph& g, edgeCountMap& triangles);
	int getOverlap(const RankedNeighbors& first, const RankedNeighbors& second);
};

} /* namespace NetworKit */
#endif /* SIMMELIANBACKBONE_H_ */
