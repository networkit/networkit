/*
 * IndependentSetFinder.h
 *
 *  Created on: 27.02.2013
 *      Author: cls
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
	virtual std::vector<bool> run(Graph& G) = 0;

	virtual std::string toString() const;
};

} /* namespace EnsembleClustering */
#endif /* INDEPENDENTSETFINDER_H_ */
