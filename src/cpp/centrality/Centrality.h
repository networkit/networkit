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
 * Abstract base class for centrality measures.
 */
class Centrality {
public:
	Centrality(const Graph& G, bool normalized=false);

	virtual ~Centrality() = default;

	virtual void run() = 0;

	virtual std::vector<double> scores();

	virtual std::vector<std::pair<node, double> > ranking();

	virtual double score(node v);

protected:

	const Graph& G;
	std::vector<double> scoreData;
	bool normalized; // true if scores should be normalized in the interval [0,1]


};

} /* namespace NetworKit */

#endif /* CENTRALITY_H_ */
