/*
* ApproxHopPlot.h
*
*  Created on: 30.03.2016
*      Author: Maximilian Vogel
*/

#ifndef APPROXHOPPLOT_H_
#define APPROXHOPPLOT_H_

#include <map>
#include "../graph/Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup distance
 */
class ApproxHopPlot : public Algorithm {

	/*
	these are variatons of the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining
	in Massive Graphs" by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
	*/
public:
	/**
	* computes the hop-plot of a given graph
	* the hop plot is the set of pairs (d, g(g)) for each natural number d and where g(d) is the fraction of connected node pairs whose shortest connecting path has length at most d
	* @param G the given graph
	* @param maxDistance the maximum path length that shall be considered. set 0 for infinite
	* @param k the number of parallel approximations to get a more robust result
	* @param r the amount of bits that are added to the length of the bitmask to improve the accuracy
	* @return the approximated hop-plot of the graph
	*/
	ApproxHopPlot(const Graph& G, const count maxDistance=0, const count k=64, const count r=7);

	void run() override;

	std::map<count, double> getHopPlot() const;

private:
	const Graph& G;
	const count maxDistance;
	const count k;
	const count r;
	std::map<count, double> hopPlot;

};

} /* namespace NetworKit */

#endif /* APPROXHOPPLOT_H_ */
