/*
 * Centrality.h
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef CENTRALITY_H_
#define CENTRALITY_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Abstract base class for centrality measures.
 */
class Centrality {
public:
	/**
	 * Constructs the Centrality class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to @c true.
	 *
	 * @param G The graph.
	 * @param normalized If set to @c true the scores are normalized in the interval [0,1].
	 */
	Centrality(const Graph& G, bool normalized=false);

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

protected:

	const Graph& G;
	std::vector<double> scoreData;
	bool normalized; // true if scores should be normalized in the interval [0,1]


};

} /* namespace NetworKit */

#endif /* CENTRALITY_H_ */
