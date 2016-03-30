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
	* computes the effective diameter exactly
	* @param G the given graph
	* @param ratio the ratio of nodes that should be connected (0,1]
	* @return the exact effective diameter of the graph
	*/
	EffectiveDiameter(const Graph& G, const double ratio=0.9);

	void run() override;

	double getEffectiveDiameter() const;

private:
	const Graph& G;
	const double ratio;
	double effectiveDiameter;
};

} /* namespace NetworKit */

#endif /* EFFECTIVEDIAMETER_H_ */
