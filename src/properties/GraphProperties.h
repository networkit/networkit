/*
 * GraphProperties.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef GRAPHPROPERTIES_H_
#define GRAPHPROPERTIES_H_


#include "../graph/Graph.h"
#include "../graph/BFS.h"
#include "../io/METISGraphReader.h"
#include "ApproximateClusteringCoefficient_Hoske.h"
#include "../graph/BFS.h"
#include "../graph/GraphDistance.h"

namespace NetworKit {

/**
 * Collection of methods for basic network properties.
 */
class GraphProperties {
public:
	GraphProperties();
	virtual ~GraphProperties();

	static std::vector<count> degreeDistribution(const Graph& G);


	/**
	 * The local clustering coefficient for a node is the number of edges among its
	 * neighbors divided by the number of possible edges.
	 * 	For an undirected graph where N(v) does not include v itself:
	 * 		$c_v := \frac{2 \cdot |E(N(v))| }{\deg(v) ( \deg(v) - 1 )}$
	 *
	 *
	 * @param[in]	G	the graph
	 * @param[out]		node -> local clustering coefficient
	 */
	static std::vector<double> localClusteringCoefficients(const Graph& G);


	/**
	 * The average local clustering coefficient for the graph.
	 * 		$\frac{1}{n} \cdot \sum_{v \in V} c_v$
	 *
	 * @param[in]	G	the graph
	 */
	static double averageLocalClusteringCoefficient(const Graph& G);

	static std::vector<double> localClusteringCoefficientPerDegree(const Graph& G);

	static std::pair<count, count> minMaxDegree(const Graph& G);

	static double averageDegree(const Graph& G);

	static double approximateGlobalClusteringCoefficient(const Graph& G, double approximationError, double probabilityError);

	/**
	 * Estimates a range for the diameter of @a G. Based on the algorithm suggested in
	 * C. Magnien, M. Latapy, M. Habib: Fast Computation of Empirically Tight Bounds for
	 * the Diameter of Massive Graphs. Journal of Experimental Algorithmics, Volume 13, Feb 2009.
	 *
	 * @return Pair of lower and upper bound for diameter.
	 */
	static std::pair<count, count> estimatedDiameterRange(const Graph& G);
	static std::pair<count, count> estimatedDiameterRange_Hoske(const Graph& G, double error);

    /**
     * Computes the ecc_G(u) and returns a vertex v with largest distance from u. 
     *
     * @return pair of node v and eccentricity.
     */
    static std::pair<node, count> eccentricity_Hoske(const Graph& G, node u);

    static std::pair<count, count> estimateDiameter_ck(const Graph& G);

    static std::pair<count, count> estimatedDiameterRange_Feist(const Graph& G, count k);
    static count DiameterRange_Feist(const Graph& G);

    static std::pair<count, count> estimatedDiameterRange_Brueckner(const Graph& G);
	static count exactDiameter_Brueckner(const Graph& G);

	static std::pair<count, count> estimatedDiameterRange_Ritter(const Graph& G);
};

} /* namespace NetworKit */
#endif /* GRAPHPROPERTIES_H_ */
