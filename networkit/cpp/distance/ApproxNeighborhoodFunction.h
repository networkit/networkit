/*
* ApproxNeighborhoodFunction.h
*
*  Created on: 30.03.2016
*      Author: Maximilian Vogel
*/

#ifndef APPROXNEIGHBORHOODFUNCTION_H_
#define APPROXNEIGHBORHOODFUNCTION_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup distance
 */
class ApproxNeighborhoodFunction : public Algorithm {

public:
	/**
	* approximates the effective diameter of a given graph
	* the effective diameter equals the number of edges on average to reach 90% of all other nodes
	*
	* @param G the given graph
	* @param ratio the ratio of nodes that should be connected (0,1]
	* @param k the number of parallel approximations to get a more robust result
	* @param r the amount of bits that are added to the length of the bitmask to improve the accuracy
	* @return the approximated effective diameter of the graph
	*/
	ApproxNeighborhoodFunction(const Graph& G, const count k=64, const count r=7);

	void run() override;

	std::vector<count> getNeighborhoodFunction() const;

private:
	const Graph& G;
	const count k;
	const count r;
	std::vector<count> result;
};

} /* namespace NetworKit */

#endif /* APPROXNEIGHBORHOODFUNCTION_H_ */
