/*
* NeighborhoodFunction.h
*
*  Created on: 30.03.2016
*      Author: Maximilian Vogel
*/

#ifndef NEIGHBORHOODFUNCTION_H_
#define NEIGHBORHOODFUNCTION_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup distance
 */
class NeighborhoodFunction : public Algorithm {

public:
	/**
	* Computes the neighborhood function exactly.
	* The neighborhood function N of a graph G for a given distance t is defined
	* as the number of node pairs (u,v) that can be reached within distance t.
	*
	* @param G the given graph
	* @return the exact effective diameter of the graph
	*/
	NeighborhoodFunction(const Graph& G);

	void run() override;

	/**
	 * Returns the neighborhood function of the graph.
	 * @return the neighborhood function of the graph
	 */
	std::vector<count> getNeighborhoodFunction() const;

private:
	const Graph& G;
	std::vector<count> result;
};

} /* namespace NetworKit */

#endif /* NEIGHBORHOODFUNCTION_H_ */
