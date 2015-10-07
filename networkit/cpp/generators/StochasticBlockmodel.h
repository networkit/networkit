/*
* StochasticBlockmodel.h
*
*  Created on: 13.08.2014
*      Author: Christian Staudt
*/

#ifndef STOCHASTICBLOCKMODEL_H_
#define STOCHASTICBLOCKMODEL_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {


class StochasticBlockmodel: public StaticGraphGenerator {

public:
	/**
	* Construct a undirected regular ring lattice.
	*
	* @param nNodes 		number of nodes in target graph
	* @param n		number of blocks (=k)
	* @param membership		maps node ids to block ids (consecutive, 0 <= i < nBlocks)
	* @param affinity		matrix of size k x k with edge probabilities betweeen the blocks
	*/
	StochasticBlockmodel(count n, count nBlocks, const std::vector<index>& membership, const std::vector<std::vector<double> >& affinity);

	virtual Graph generate();

protected:
		count n;
		count nBlocks;
		std::vector<index> membership;
		std::vector<std::vector<double> > affinity;

};

} /* namespace NetworKit */
#endif /* STOCHASTICBLOCKMODEL_H_ */
