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
	* computes the effective diameter exactly
	* @param G the given graph
	* @param ratio the ratio of nodes that should be connected (0,1]
	* @return the exact effective diameter of the graph
	*/
	NeighborhoodFunction(const Graph& G);

	void run() override;

	std::vector<count> getNeighborhoodFunction() const;

private:
	const Graph& G;
	std::vector<count> result;
};

} /* namespace NetworKit */

#endif /* NEIGHBORHOODFUNCTION_H_ */
