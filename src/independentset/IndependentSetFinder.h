/*
 * IndependentSetFinder.h
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef INDEPENDENTSETFINDER_H_
#define INDEPENDENTSETFINDER_H_

#include "../graph/Graph.h"

namespace EnsembleClustering {


class IndependentSetFinder {


public:

	IndependentSetFinder();

	virtual ~IndependentSetFinder();

	/**
	 * @return a boolean vector of length n where vec[v] is true iff v is in the independent sets
	 *
	 * @param[in]	G	graph
	 */
	virtual std::vector<bool> run(const Graph& G) = 0;

	virtual std::string toString() const;

	/**
	 * Check whether a set is independent.
	 */
	bool isIndependentSet(const std::vector<bool>& set, const Graph& G) const;

};

} /* namespace EnsembleClustering */
#endif /* INDEPENDENTSETFINDER_H_ */
