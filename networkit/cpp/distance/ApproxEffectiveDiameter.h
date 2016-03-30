/*
* ApproxEffectiveDiameter.h
*
*  Created on: 29.03.16
*      Author: Maximilian Vogel
*/

#ifndef APPROXEFFECTIVEDIAMETER_H_
#define APPROXEFFECTIVEDIAMETER_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup distance
 */
class ApproxEffectiveDiameter : public Algorithm {

	/*
	these are variatons of the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining
	in Massive Graphs" by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
	*/
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
	ApproxEffectiveDiameter(const Graph& G, const double ratio=0.9, const count k=64, const count r=7);

	void run() override;

	double getEffectiveDiameter() const;

private:
	const Graph& G;
	const double ratio;
	const count k;
	const count r;
	double effectiveDiameter;
};

} /* namespace NetworKit */

#endif /* APPROXEFFECTIVEDIAMETER_H_ */
