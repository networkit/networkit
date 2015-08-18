/*
 * EdgeScore.h
 *
 *  Created on: 18.08.2015
 *      Author: Gerd Lindner
 */

#ifndef EDGESCORE_H_
#define EDGESCORE_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"
#include <vector>

namespace NetworKit {
/**
 * Abstract base class for an edge score.
 */
template<typename T>
class EdgeScore  : public Algorithm {

public:

	EdgeScore(const Graph& G) : Algorithm(), G(G) {
		if (G.isDirected()) {
			WARN("Application to directed graphs is not well tested");
		}
	}

	/** Compute the edge score. */
	virtual void run() {
		// empty run method for edge scoring methods that do not require preprocessing but calculate score(u,v) on the fly
		hasRun = true;
	};

	/** Get a vector containing the score for each edge in the graph.
	@Return the edge scores calculated by @link run().
	*/
	virtual std::vector<T> scores() const {
		if (!hasRun) {
			throw std::runtime_error("Call run method first");
		}
		return scoreData;
	}

	/** Get the edge score of the edge with the given edge id.
	*/
	virtual T score(edgeid eid) = 0;

	/** Get the edge score of the given edge.
	*/
	virtual T score(node u, node v) = 0;

protected:
	const Graph& G;
	std::vector<T> scoreData;


};

}


#endif /* EDGESCORE_H_ */
