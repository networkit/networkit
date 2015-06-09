/*
 * Centrality.h
 *
 *  Created on: 19.02.2014
 *      Author: Christian Staudt
 */

#ifndef CENTRALITY_H_
#define CENTRALITY_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Abstract base class for centrality measures.
 */
class Centrality : public Algorithm {
public:
	/**
	 * Constructs the Centrality class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to @c true.
	 *
	 * @param G The graph.
	 * @param normalized If set to @c true the scores are normalized in the interval [0,1].
	 * @param computeEdgeCentrality		If true, compute also edge centralities (for algorithms where this is applicable)
	 */
	Centrality(const Graph& G, bool normalized=false, bool computeEdgeCentrality=false);

	/** Default destructor */
	virtual ~Centrality() = default;

	/**
	 * Compute betweenness scores.
	 */
	virtual void run() = 0;

	/**
	 * Get a vector containing the betweenness score for each node in the graph.
	 * @return The betweenness scores calculated by @link run().
	 */
	virtual std::vector<double> scores();

	/**
	 * Get a vector containing the edge betweenness score for each edge in the graph.
	 * @return The edge betweenness scores calculated by @link run().
	 */
	virtual std::vector<double> edgeScores();

	/**
	 * Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score
	 * calculated by @link run().
	 * @return A vector of pairs.
	 */
	virtual std::vector<std::pair<node, double> > ranking();

	/**
	 * Get the betweenness score of node @a v calculated by @link run().
	 *
	 * @param v A node.
	 * @return The betweenness score of node @a v.
	 */
	virtual double score(node v);

	/**
	* Get the theoretical maximum of centrality score in the given graph.
	*
	* @return The maximum centrality score.
	*/
	virtual double maximum();

protected:

	const Graph& G;
	std::vector<double> scoreData;
	std::vector<double> edgeScoreData;
	bool normalized; // true if scores should be normalized in the interval [0,1]
	bool computeEdgeCentrality;

};

} /* namespace NetworKit */

#endif /* CENTRALITY_H_ */
