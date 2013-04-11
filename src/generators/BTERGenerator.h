/*
 * BTERGenerator.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef BTERGENERATOR_H_
#define BTERGENERATOR_H_

#include <algorithm>
#include <cmath>

#include "StaticGraphGenerator.h"

namespace NetworKit {

typedef count degree; //!< node degree is the number of incident edges on a node

/**
 * Block Two-level Erdos-Renyi Graph Model
 *
 * Paper: http://www.sandia.gov/~tgkolda/pubs/bibtgkfiles/BTER-arXiv-1302.6636v1.pdf
 * Video introduction: http://www.youtube.com/watch?v=kF-pKfg846M
 *
 *
 */
class BTERGenerator: public NetworKit::StaticGraphGenerator {


/**
 * Definitions
 *
 * affinity block:
 * 		nodes within the same affinity block have a much higher chance of being connected than nodes at random.
 * affinity group:
 * 		all affinity blocks of the same size and minimum degree can be grouped together into an affinity group Ñ
 * 		all blocks in the same group share the same block size and weight.
 * bulk node:
 * 		In a block where all nodes are the same degree, we say the nodes are bulk nodes.
 * filler node:
 * 		In a block with nodes of differing degrees, all nodes with degree equal to the minimum degree are still bulk nodes.
 * 		The remaining nodes are called filler nodes.
 *
 */



public:

	/**
	 * @param[in]	degreeDistribution			number of nodes per degree
	 * @param[in]	clusteringCoefficients		mean clusering coefficient for nodes per degree
	 */
	BTERGenerator(std::vector<count> degreeDistribution, std::vector<count> clusteringCoefficients);

	virtual ~BTERGenerator();

	virtual Graph generate();


protected:

	void preprocessing();

	void phaseOne();

	void phaseTwo();

	void BTERSample();

	void BTERSamplePhaseOne();

	void BTERSamplePhaseTwo();

	void BTERSamplePhase2Node();


	std::vector<count> n_; //!< degree distribution: n_[d] = number of nodes with degree d
	std::vector<double> c_;	//!< clustering coefficient per degree:
	// TODO: no need to keep c_ after the preprocessing is done
	count dMax; //!< maximum degree

};


} /* namespace NetworKit */
#endif /* BTERGENERATOR_H_ */
