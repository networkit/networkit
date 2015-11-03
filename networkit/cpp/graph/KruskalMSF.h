/*
 * KruskalMSF.h
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#ifndef KRUSKALMSF_H_
#define KRUSKALMSF_H_

#include "../Globals.h"
#include "Graph.h"
#include "SpanningForest.h"

namespace NetworKit {

/**
 * Creates a minimum spanning tree for each connected component.
 * @ingroup graph
 */
class KruskalMSF: public SpanningForest {
public:
	KruskalMSF(const Graph& G);
	virtual ~KruskalMSF() = default;

	/**
	 * Computes for each component a minimum weight spanning tree
	 * (or simply a spanning tree in unweighted graphs).
	 * Uses Kruskal's algorithm.
	 * Time complexity: sort(n) + n * inverse Ackermann(n, m).
	 */
	virtual void run() override;
};

} /* namespace NetworKit */
#endif /* KRUSKALMSF_H_ */
