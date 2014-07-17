/*
 * Diameter.h
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#ifndef DIAMETER_H_
#define DIAMETER_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup properties
 */
class Diameter {

public:

	/**
	 * Estimates a range for the diameter of @a G. Based on the algorithm suggested in
	 * C. Magnien, M. Latapy, M. Habib: Fast Computation of Empirically Tight Bounds for
	 * the Diameter of Massive Graphs. Journal of Experimental Algorithmics, Volume 13, Feb 2009.
	 *
	 * @return Pair of lower and upper bound for diameter.
	 */
	static std::pair<edgeweight, edgeweight> estimatedDiameterRange(const Graph& G, double error);

	/**
	 * Get the exact diameter of the graph @a G. The algorithm for unweighted graphs is based on
	 * the iFub algorithm suggested in Pilu Crescenzi, Roberto Grossi, Michel Habib, Leonardo Lanzi, Andrea Marino:
	 * On computing the diameter of real-world undirected graphs,
	 * Theoretical Computer Science, Volume 514, 25 November 2013, Pages 84-95, ISSN 0304-3975,
	 * http://dx.doi.org/10.1016/j.tcs.2012.09.018
	 *
	 * @param G The graph.
	 * @return exact diameter of the graph @a G
	 */
	static edgeweight exactDiameter(const Graph& G);


	/**
	 * Get a 2-approximation of the node diameter (unweighted diameter) of @a G.
	 *
	 * @param[in]	G			The graph.
	 * @param[in]	samples		One sample is enough if the graph is connected. If there
	 *							are multiple connected components, then the number of samples
	 *							must be chosen so that the probability of sampling the component
	 *							with the largest diameter ist high.
	 * @return A 2-approximation of the vertex diameter (unweighted diameter) of @a G.
	 */
	static edgeweight estimatedVertexDiameter(const Graph& G, count samples);


	/** @return a 2-approximation of the vertex diameter (unweighted diameter) of @a G.
			Considers each connected component and returns the maximum diameter.
	 */
	static edgeweight estimatedVertexDiameterPedantic(const Graph& G);
};

} /* namespace NetworKit */

#endif /* DIAMETER_H_ */
