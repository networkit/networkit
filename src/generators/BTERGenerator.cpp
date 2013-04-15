/*
 * BTERGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "BTERGenerator.h"

namespace NetworKit {

BTERGenerator::BTERGenerator(std::vector<count>& degreeDistribution,
		std::vector<double>& clusteringCoefficients, double beta) :
		beta(beta),
		dMax(degreeDistribution.size() - 1), // degree distribution has entries for indices in [0, dMax]
		nd_(degreeDistribution),
		c_(clusteringCoefficients),
		id_(dMax),
		nFill_(dMax),
		wd_(dMax),
		r_(dMax, none),
		ig_(dMax, none), 	// gMax <= dMax TODO: save memory here?
		b_(dMax, none),	// gMax <= dMax
		wg_(dMax, none),	// gMax <= dMax
		ng_(dMax, none) 	// gMax <= dMax
{
}

BTERGenerator::~BTERGenerator() {
	// TODO Auto-generated destructor stub
}

std::pair<count, count> BTERGenerator::desiredGraphSize() {
	assert (nd_.size() == dMax + 1);
	count nDesired = 0;
	count mDesired = 0;
	for (index d = 0; d <= dMax; ++d) {
		nDesired += nd_[d];
		mDesired += d * nd_[d];
	}
	mDesired = mDesired / 2;
	return std::pair<count, count>(nDesired, mDesired);
}

Graph BTERGenerator::generate() {

	std::pair<count, count> size = this->desiredGraphSize();
	INFO("desired number of nodes: " << size.first);
	INFO("desired number of edges: " << size.second);

	DEBUG("setup");
	this->setup();
	DEBUG("sample");
	this->sample(); // TODO: insert edges directly



	return Graph(0); // FIXME: mock graph
}

void BTERGenerator::setup() {

	// assign nodes to affinity blocks of d +1 nodes with degree d

	// compute number of nodes n'_d with degree greater d

	std::vector<count> nBulk_(dMax); // number of bulk nodes of degree d
	std::vector<count> ndRest_(dMax); // ?
	std::vector<double> wFill_(dMax); // weight of the degree-d fill nodes for phase 2
	std::vector<double> wBulk_(dMax); //weight of degree-d bulk nodes for phase 2

	// number of nodes from least degree to greatest, except degree-1 nodes are last
	id_[2] = 1; // TODO: should indices be zero-based?
	for (count d = 3; d <= dMax; d++) {
		id_[d] = id_[d - 1] + nd_[d - 1];
	}
	id_[1] = id_[dMax] + nd_[dMax];

	// compute number of nodes with degree greater than d
	for (count d = 1; d < dMax; d++) {
		for (count d_ = d + 1; d_ <= dMax; d_++) {
			ndRest_[d] += nd_[d_]; // TODO: check that higherDegreeNodes[dMax]Ê== 0
		}
	}

	// handle degree-1 nodes
	nFill_[1] = beta * nd_[1];
	wd_[1] = 0.5 * nFill_[1];
	r_[1] = 1; // TODO: should rFill_ be a seperate array from r_?

	// main loop
	count g = -1; // affinity group index - first group has index (correct?)
	count nFillStar = 0; // number of nodes needed to complete the last incomplete block
	count dStar = 0; // internal degree of the last incomplete block
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
			ig_[g] = id_[d] + nFill_[d];
			b_[g] = std::ceil(nBulk_[d] / (d + 1));
			nd_[g] = d + 1;
			if ((b_[g] * (d + 1)) > (ndRest_[d] + nBulk_[d])) {
				// special handling of last group
				if (b_[g] != 1) {
					throw std::runtime_error("?"); // TODO: what happens here? check for indexing error
				}
				nd_[g] = (ndRest_[d] + nBulk_[d]);
			}
			rhoStar = std::pow(c_[d], (1.0 / 3.0)); //
			dStar = (nd_[g] - 1) * rhoStar;
			wBulk_[d] = 0.5 * nBulk_[d] * (d - dStar);
			DEBUG("wg_: " << wg_.size() << " b_: " << b_.size() << " ng_: " << ng_.size());
			assert ( (1 - rhoStar) > 0); // avoid division by 0
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
	double w1 = 0.0;
	for (double wg : wg_) {
		w1 += wg;
	}
	double w2 = 0.0;
	for (double wd : wd_) {
		w2 += wd;
	}
	double w = w1 + w2;

	std::vector<std::pair<node, node> > E1;
	std::vector<std::pair<node, node> > E2;

	DEBUG("start sampling");
	for (count j = 0; j < w; ++j) {
		double r = this->rand.probability();
		if (r < (w1 / w)) {
			E1.push_back(this->samplePhaseOne());
		} else {
			E2.push_back(this->samplePhaseTwo());
		}
	}

	INFO("E_1 size: " << E1.size());
	INFO("E_2 size: " << E2.size());
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
