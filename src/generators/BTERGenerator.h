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
#include "../aux/Random.h"

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


public:

	/**
	 * @param[in]	degreeDistribution			number of nodes per degree
	 * @param[in]	clusteringCoefficients		mean clusering coefficient for nodes per degree
	 * @param[in]	beta						blowup factor for degree-1 nodes
	 */
	BTERGenerator(std::vector<count>& degreeDistribution, std::vector<double>& clusteringCoefficients, double beta = 1.0);

	virtual ~BTERGenerator();

	virtual Graph generate();


protected:

	void setup();

	void sample();

	std::pair<node, node> samplePhaseOne();

	std::pair<node, node> samplePhaseTwo();

	node samplePhaseTwoNode();



	Aux::Random rand; // random module

	double beta; // blowup-factor for deg-1 nodes

	degree dMax; //!< maximum degree
	std::vector<count> nd_; //!< degree distribution: nd_[d] = number of nodes with degree d
	std::vector<double> c_;	//!< clustering coefficient per degree: 	// TODO: no need to keep c_ after the preprocessing is done

	// instead of passing these as arguments, keep them as member variables

	std::vector<index> id_; // index i_d for first node of each degree d
	std::vector<count> nFill_; // ?
	std::vector<double> wd_; // ? double or int? sum of wFill and wBulk
	std::vector<double> r_; // ratio of fill nodes
	std::vector<index> ig_; //? group index vector: ?
	std::vector<count> b_; // b_[g]: number of blocks in a group
	std::vector<double> wg_; // wg_[g]: weight of the group g
	std::vector<count> ng_; // ?






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


};


} /* namespace NetworKit */
#endif /* BTERGENERATOR_H_ */
