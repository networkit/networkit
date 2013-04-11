/*
 * BTERGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "BTERGenerator.h"

namespace NetworKit {

BTERGenerator::BTERGenerator(std::vector<count> degreeDistribution,
		std::vector<count> clusteringCoefficients) : n_(degreeDistribution), dMax(n_.size()) {
}

BTERGenerator::~BTERGenerator() {
	// TODO Auto-generated destructor stub
}

Graph BTERGenerator::generate() {
}

void BTERGenerator::preprocessing() {

	double beta; // blowup-factor for deg-1 nodes TODO: parameter

	// TODO: assign nodes to affinity blocks of d +1 nodes with degree d

	// TODO: compute number of nodes n'_d with degree greater d

	// number of nodes from least degree to greatest, except degree-1 nodes are last
	std::vector<index> id_; // index i_d for first node of each degree d
	id_[2] = 1; // TODO: should indices be zero-based?
	for (count d = 3; d <= dMax; d++) {
		id_[d] = id_[d - 1] + n_[d - 1];
	}
	id_[1] = id_[dMax] + n_[dMax];

	// compute number of nodes with degree greater than d
	std::vector<count> higherDegreeNodes;
	for (count d = 1; d < dMax; d++) {
		for (count d_ = d + 1; d_ <= dMax; d_++) {
			higherDegreeNodes[d] += n_[d_]; // TODO: check that higherDegreeNodes[dMax]Ê== 0
		}
	}

	// handle degree-1 nodes
	count n1Fill = beta * n_[1];
	double w1 = 0.5 * n1Fill;
	count r1Fill = 1;

	// main loop
	count g = 0; // affinity group index ?
	count nFillStar = 0; // ?
	count dStar = 0; // ?
	std::vector<count> nFill_; // ?
	std::vector<double> wFill_; // ?
	std::vector<count> nBulk_; // ?
	std::vector<count> b_; // ?  index vector?
	std::vector<count> nRest_; // ?
	std::vector<double> wBulk_; // ?
	double rhoStar; // ?

	std::vector<index> ig_; //? group index vector: ?

	std::vector<count> ng_; // ?
	std::vector<double> w_; // group weight ?

	std::vector<double> r_; // ?


	for (count d = 2; d <= dMax; d++) {
		// try to fill incomplete block from current group
		if (nFillStar > 0) {
			nFill_[d] = std::min(nFillStar, n_[d]);
			nFillStar = nFillStar - nFill_[d];
			wFill_[d] = 0.5 * nFill_[d] * (d - dStar);
		} else {
			nFill_[d] = 0;
			wFill_[d] = 0.0;
		}
		nBulk_[d] = n_[d] - nFill_[d];

		if (nBulk_[d] > 0) {
			// create a new group for degree-d bulk nodes
			g += 1;
			ig_[g] = id_[d] + nFill_[d]; // ? is i_g and i_d the same array?
			b_[g] = std::ceil(nBulk_[d] / (d + 1));
			n_[g] = d + 1;
			if ((b_[g] * (d + 1)) > (nRest_[d] + nBulk_[d])) {
				// special handling of last group
				if (b_[g] != 1) {
					throw std::runtime_error("?");
				}
				n_[g] = (nRest_[d] + nBulk_[d]);
			}
			rhoStar = std::pow(c_[d], (1.0 / 3.0)); //
			dStar = (n_[g] - 1) * rhoStar;
			wBulk_[d] = 0.5 * nBulk_[d] * (d - dStar);
			w_[g] = b_[g] * 0.5 * ng_[g] * (ng_[g] - 1) * std::log(1.0 / (1 - rhoStar)); // correct log?
			nFillStar = (b_[g] * ng_[g]) - nBulk_[d];
		} else {
			wBulk_[d] = 0;
		}
		w_[d] = wFill_[d] + wBulk_[d];
		r_[d] = wFill_[d] / w_[d];
	}

	// arrays are passed to other procesures as member variables
}

void BTERGenerator::phaseOne() {
	// TODO: requires: mapping of nodes to affinity blocks - do intelligent implementation
}

void BTERGenerator::phaseTwo() {
	// TODO: requires: expected excess degree of
}

void BTERGenerator::BTERSample() {
}

void BTERGenerator::BTERSamplePhaseOne() {
}

void BTERGenerator::BTERSamplePhaseTwo() {
}

void BTERGenerator::BTERSamplePhase2Node() {
	// degree d = w_[];
}

} /* namespace NetworKit */
