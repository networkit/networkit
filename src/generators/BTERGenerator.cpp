/*
 * BTERGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "BTERGenerator.h"

namespace NetworKit {

BTERGenerator::BTERGenerator(std::vector<count> degreeDistribution,
		std::vector<count> clusteringCoefficients) : nd_(degreeDistribution), dMax(nd_.size()) {
}

BTERGenerator::~BTERGenerator() {
	// TODO Auto-generated destructor stub
}

Graph BTERGenerator::generate() {
	this->setup();
}

void BTERGenerator::setup() {


	std::vector<double> wFill_; // ?
	std::vector<count> nBulk_; // ?
	std::vector<count> nRest_; // ?
	std::vector<double> wBulk_; // ?
	std::vector<double> rFill_;
	std::vector<double> r_; // ?

	// TODO: assign nodes to affinity blocks of d +1 nodes with degree d

	// TODO: compute number of nodes n'_d with degree greater d

	// number of nodes from least degree to greatest, except degree-1 nodes are last
	id_[2] = 1; // TODO: should indices be zero-based?
	for (count d = 3; d <= dMax; d++) {
		id_[d] = id_[d - 1] + nd_[d - 1];
	}
	id_[1] = id_[dMax] + nd_[dMax];

	// compute number of nodes with degree greater than d
	std::vector<count> higherDegreeNodes;
	for (count d = 1; d < dMax; d++) {
		for (count d_ = d + 1; d_ <= dMax; d_++) {
			higherDegreeNodes[d] += nd_[d_]; // TODO: check that higherDegreeNodes[dMax]Ê== 0
		}
	}

	// handle degree-1 nodes
	nFill_[1] = beta * nd_[1];
	wd_[1] = 0.5 * nFill_[1];
	rFill_[1] = 1;

	// main loop
	count g = 0; // affinity group index ?
	count nFillStar = 0; // ?
	count dStar = 0; // ?
	double rhoStar; // ?

	for (count d = 2; d <= dMax; d++) {
		// try to fill incomplete block from current group
		if (nFillStar > 0) {
			nFill_[d] = std::min(nFillStar, nd_[d]);
			nFillStar = nFillStar - nFill_[d];
			wFill_[d] = 0.5 * nFill_[d] * (d - dStar);
		} else {
			nFill_[d] = 0;
			wFill_[d] = 0.0;
		}
		nBulk_[d] = nd_[d] - nFill_[d];

		if (nBulk_[d] > 0) {
			// create a new group for degree-d bulk nodes
			g += 1;
			ig_[g] = id_[d] + nFill_[d]; // ? is i_g and i_d the same array?
			b_[g] = std::ceil(nBulk_[d] / (d + 1));
			nd_[g] = d + 1;
			if ((b_[g] * (d + 1)) > (nRest_[d] + nBulk_[d])) {
				// special handling of last group
				if (b_[g] != 1) {
					throw std::runtime_error("?");
				}
				nd_[g] = (nRest_[d] + nBulk_[d]);
			}
			rhoStar = std::pow(c_[d], (1.0 / 3.0)); //
			dStar = (nd_[g] - 1) * rhoStar;
			wBulk_[d] = 0.5 * nBulk_[d] * (d - dStar);
			wg_[g] = b_[g] * 0.5 * ng_[g] * (ng_[g] - 1) * std::log(1.0 / (1 - rhoStar)); // correct log?
			nFillStar = (b_[g] * ng_[g]) - nBulk_[d];
		} else {
			wBulk_[d] = 0;
		}
		wd_[d] = wFill_[d] + wBulk_[d];
		r_[d] = wFill_[d] / wd_[d];
	}

	// arrays are passed to other procesures as member variables
}


void BTERGenerator::sample() {
}

std::pair<node, node> BTERGenerator::samplePhaseOne() {
	index g = this->rand.choice(wd_); // choose group
	double r1 = this->rand.probability();
	double r2 = this->rand.probability();
	double r3 = this->rand.probability();
	index delta = ig_[g] + std::floor(r1 * b_[g]) * ng_[g]; // choose block and compute its offset
	node u = std::floor(r2 * ng_[g]) + delta;
	node v = std::floor(r3 * (ng_[g] - 1)) + delta;
	if (v >= u) {
		v += 1; // TODO: why?
	}
	return std::pair<node, node>(u, v); // TODO: create edge here?
}

std::pair<node, node> BTERGenerator::samplePhaseTwo() {
	node u = this->samplePhaseTwoNode();
	node v = this->samplePhaseTwoNode();
	return std::pair<node, node>(u, v);
}

node BTERGenerator::samplePhaseTwoNode() {
	degree d = this->rand.choice(wd_);
	double r1 = this->rand.probability();
	double r2 = this->rand.probability();
	node u;
	if (r1 < nFill_[d]) {
		u = std::floor(r2 * nFill_[d]) + id_[d];
	} else {
		u = std::floor(r2 * (nd_[d] - nFill_[d])) + (id_[d] + nFill_[d]);
	}
	return u;
}

} /* namespace NetworKit */
