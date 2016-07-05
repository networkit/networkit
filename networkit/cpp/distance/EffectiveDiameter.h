/*
* EffectiveDiameter.h
*
*  Created on: 16.06.2014
*      Author: Marc Nemes
*/

#ifndef EFFECTIVEDIAMETER_H_
#define EFFECTIVEDIAMETER_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup distance
 */
class EffectiveDiameter : public Algorithm {

public:
	/**
	* Computes the effective diameter exactly.
	* The effective diameter is defined as the number of edges on average to reach \p ratio of all other nodes.
	* @param G the given graph
	* @param ratio the ratio of nodes that should be connected (0,1], default = 0.9
	*/
	EffectiveDiameter(const Graph& G, const double ratio=0.9);

	void run() override;

	/**
	 * Returns the exact effective diameter of the graph.
	 * @return the exact effective diameter of the graph
	 */
	double getEffectiveDiameter() const;

private:
	const Graph& G;
	const double ratio;
	double effectiveDiameter;
};

} /* namespace NetworKit */

#endif /* EFFECTIVEDIAMETER_H_ */
