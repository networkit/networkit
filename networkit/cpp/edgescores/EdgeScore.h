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
#include <unordered_map>
#include <vector>

namespace NetworKit {
/**
 * Abstract base class for an edge score.
 */
template<typename T>
class EdgeScore  : public Algorithm {

public:

	virtual ~EdgeScore() = default;

	/** Compute the edge score. */
	virtual void run() = 0;

	/** Get a vector containing the score for each edge in the graph.
	@Return the edge scores calculated by @link run().
	*/
	virtual std::vector<T> getEdgeScores() = 0;

	std::vector<T>* _getEdgeScores() {
		return new std::vector<T>{std::move(getEdgeScores())};
	};

};

}


#endif /* EDGESCORE_H_ */
