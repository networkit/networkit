/*
* NeighborhoodFunctionHeuristic.h
*
*      Author: Maximilian Vogel
*/

#ifndef NEIGHBORHOODFUNCTIONHEURISTIC_H_
#define NEIGHBORHOODFUNCTIONHEURISTIC_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup distance
 */
class NeighborhoodFunctionHeuristic : public Algorithm {

public:
	enum SelectionStrategy {
		RANDOM,
		SPLIT
	};
	/**
	* Computes a heuristic of the neighborhood function.
	*
	* The algorithm runs nSamples breadth-first searches and scales the results up to the actual amount of nodes.
	* Accepted strategies are "split" and "random".
	*
	* @param G the given graph
	* @param nSamples the amount of samples, set to zero for heuristic of max(sqrt(m), 0.15*n)
	* @param strategy the strategy to select the samples, accepts "random" or "split"
	*/
	NeighborhoodFunctionHeuristic(const Graph& G, const count nSamples = 0, const SelectionStrategy strategy = SPLIT);

	void run() override;

	/**
	 * Returns the approximated neighborhood function of the graph.
	 * @return the approximated neighborhood function of the graph
	 */
	std::vector<count> getNeighborhoodFunction() const;

private:
	const Graph& G;
	const count nSamples;
	const SelectionStrategy strategy;
	std::vector<count> result;

	/* selection schemes implemented as private functions */
	std::vector<node> random(const Graph& G, count nSamples);
	std::vector<node> split(const Graph& G, count nSamples);

};

} /* namespace NetworKit */


#endif /* NEIGHBORHOODFUNCTIONHEURISTIC_H_ */