/*
* EffectiveDiameter.h
*
*  Created on: 16.06.2014
*      Author: Marc Nemes
*/

#ifndef EFFECTIVEDIAMETER_H_
#define EFFECTIVEDIAMETER_H_

#include <map>
#include "../graph/Graph.h"

namespace NetworKit {

class EffectiveDiameter {

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
	static double effectiveDiameter(const Graph& G, const double ratio=0.9, const count k=64, const count r=7);

	/**
	* computes the effective diameter exactly
	* @param G the given graph
	* @param ratio the ratio of nodes that should be connected (0,1]
	* @return the exact effective diameter of the graph
	*/
	static double effectiveDiameterExact(const Graph& G, const double ratio=0.9);

	/**
	* computes the hop-plot of a given graph
	* the hop plot is the set of pairs (d, g(g)) for each natural number d and where g(d) is the fraction of connected node pairs whose shortest connecting path has length at most d
	* @param G the given graph
	* @param maxDistance the maximum path length that shall be considered. set 0 for infinite
	* @param k the number of parallel approximations to get a more robust result
	* @param r the amount of bits that are added to the length of the bitmask to improve the accuracy
	* @return the approximated hop-plot of the graph
	*/
	static std::map<count, double> hopPlot(const Graph& G, const count maxDistance=0, const count k=64, const count r=7);

};

} /* namespace NetworKit */

#endif /* EFFECTIVEDIAMETER_H_ */
