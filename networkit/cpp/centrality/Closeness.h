/*
 * Closeness.h
 *
 *  Created on: 03.10.2014
 *      Author: nemes
 */

#ifndef CLOSENESS_H_
#define CLOSENESS_H_

#include "Centrality.h"

namespace NetworKit {

enum ClosenessVariant { standard = 0, generalized = 1};

/**
 * @ingroup centrality
 */
class Closeness : public Centrality {

public:
	/**
	 * Constructs the Closeness class for the given Graph @a G. If the closeness
	 * scores should be normalized, then set @a normalized to <code>true</code>.
	 * The run() method takes O(nm) time, where n is the number of nodes and m is
	 * the number of edges of the graph. NOTICE: the graph has to be connected.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>false</code> if scores should
	 * not be normalized into an interval of [0, 1]. Normalization only for
	 * unweighted graphs.
	 *
	 */
	Closeness(const Graph &G, bool normalized,
	          const ClosenessVariant variant = ClosenessVariant::standard);

	/**
	 * Old constructor, we keep it for backward compatibility. It computes the
	 * standard variant of the closenes.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>false</code> if scores should
	 * not be normalized into an interval of [0, 1]. Normalization only for
	 * unweighted graphs.
	 * @param checkConnectedness turn this off if you know the graph is connected.
	 *
	 */
	Closeness(const Graph &G, bool normalized = true,
	          bool checkConnectedness = true);

	/**
	 * Computes closeness cetrality on the graph passed in constructor.
	 */
	void run() override;

	/*
	 * Returns the maximum possible Closeness a node can have in a graph with the
	 * same amount of nodes (=a star)
	 */
	double maximum() override;

private:
	const count n;
	ClosenessVariant variant;
	std::vector<count> reachableNodes;

	void checkConnectedComponents() const;
};

} /* namespace NetworKit */

#endif /* CLOSENESS_H_ */
