/*
 * Contracter.h
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

	std::vector<node> getNodeMapping() const;

	virtual std::string toString() const;

protected:
	const Graph& G;
	Graph Gcoarsed;
	std::vector<node> nodeMapping;

};

} // namespace


#endif /* GRAPHCOARSENING_H_ */
