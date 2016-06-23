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
class EdgeScore : public Algorithm {

public:

	EdgeScore(const Graph& G);

	/** Compute the edge score. */
	virtual void run();

	/** Get a vector containing the score for each edge in the graph.
	@Return the edge scores calculated by @link run().
	*/
	virtual std::vector<T> scores() const;

	/** Get the edge score of the edge with the given edge id.
	*/
	virtual T score(edgeid eid);

	/** Get the edge score of the given edge.
	*/
	virtual T score(node u, node v);

protected:
	const Graph& G;
	std::vector<T> scoreData;

};

}


#endif /* EDGESCORE_H_ */
