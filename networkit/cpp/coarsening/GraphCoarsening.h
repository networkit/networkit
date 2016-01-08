/*
 * GraphCoarsening.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHCOARSENING_H_
#define GRAPHCOARSENING_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup coarsening
 * Abstract base class for graph coarsening/contraction algorithms.
 */
class GraphCoarsening : public Algorithm {

public:

	GraphCoarsening(const Graph& G);

	virtual ~GraphCoarsening() = default;

	virtual void run() = 0;

	Graph getCoarseGraph() const;

	/**
	 * Get mapping from fine to coarse node.
	 */
	std::vector<node> getFineToCoarseNodeMapping() const;

	/**
	 * Get mapping from coarse node to collection of fine nodes.
	 */
	std::map<node, std::vector<node> > getCoarseToFineNodeMapping() const;

protected:
	const Graph& G;
	Graph Gcoarsened;
	std::vector<node> nodeMapping;

};

} // namespace


#endif /* GRAPHCOARSENING_H_ */
