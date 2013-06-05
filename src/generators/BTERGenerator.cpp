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
		// array sizes are dMax + 1, indexing begins with 1
		id_(dMax + 1, none),
		nFill_(dMax + 1, none),
		wd_(dMax + 1, none),
		r_(dMax + 1, none),
		ig_(dMax + 1, none), 	// gMax <= dMax TODO: save memory here?
		b_(dMax + 1, none),	// gMax <= dMax
		wg_(dMax + 1, none),	// gMax <= dMax
		ng_(dMax + 1, none) 	// gMax <= dMax
{
	if (degreeDistribution.size() != clusteringCoefficients.size()) {
		throw std::runtime_error("degree distribution and clustering coefficients must have the same size (maximum degree + 1)");
	}
	if (beta < 1.0) {
		throw std::runtime_error("blowup factor beta must be >= 1.0");
	}

	this->G = NULL; // G is created in BTERGenerator.generate
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

	this->G = new Graph(size.first);

	INFO("BTERGenreator setup phase");
	this->setup();
	INFO("BTERGenerator sample phase");
	this->sample(); // TODO: insert edges directly

	return *(this->G);
}

void BTERGenerator::setup() {

	DEBUG("nd_ degree distribution: nd_[d] = number of nodes with degree d : " << Aux::vectorToString(nd_));
	DEBUG("c_ clustering coefficient per degree: " << Aux::vectorToString(c_));
	DEBUG("id_ : index i_d for first node of each degree d" << Aux::vectorToString(id_));
	DEBUG("nFill_  number of filler nodes per degree " << Aux::vectorToString(nFill_));
	DEBUG("wd_ sum of wFill and wBulk: " << Aux::vectorToString(wd_));
	DEBUG("r_ ratio of fill excess degree for degree : " << Aux::vectorToString(r_));
	DEBUG("ig_ start index for affinity group : " << Aux::vectorToString(ig_));
	DEBUG("b_ b_[g]: number of blocks in a group : " << Aux::vectorToString(b_));
	DEBUG("wg_ wg_[g]: weight of the group g : " << Aux::vectorToString(wg_));
	DEBUG("ng_  ng_[g]: number of blocks in the affinity group g : " << Aux::vectorToString(ng_));

	// assign nodes to affinity blocks of d +1 nodes with degree d

	// compute number of nodes n'_d with degree greater d

	std::vector<count> nBulk_(dMax + 1); 	// number of bulk nodes of degree d
	std::vector<count> ndRest_(dMax + 1); 	// ?
	// TODO: check the formula. wFill and wBulk are not not explicitly stored according to the table
	std::vector<double> wFill_(dMax + 1); 	// weight of the degree-d fill nodes for phase 2
	std::vector<double> wBulk_(dMax + 1); 	// weight of degree-d bulk nodes for phase 2

	// number of nodes from least degree to greatest, except degree-1 nodes are last
	id_[2] = 1; // TODO: should indices be zero-based?
	for (count d = 3; d <= dMax; d++) {
		id_[d] = id_[d - 1] + nd_[d - 1]; // FIXME: invalid write
	}
	id_[1] = id_[dMax] + nd_[dMax]; // FIXME: invalid read

	// compute number of nodes with degree greater than d
	for (count d = 1; d < dMax; d++) {
		for (count d_ = d + 1; d_ <= dMax; d_++) {
			ndRest_[d] += nd_[d_]; // TODO: check that higherDegreeNodes[dMax]Ê== 0
		}
	}

	// handle degree-1 nodes
	nFill_[1] = beta * nd_[1];
	wd_[1] = 0.5 * nFill_[1];
	r_[1] = 1;

	// main loop
	count g = 0; 			// affinity group index - first group has index 1 TODO: correct?
	count nFillStar = 0; 	// number of nodes needed to complete the last incomplete block
	count dStar = 0; 		// internal degree of the last incomplete block
	double rhoStar = -1; 		// ?

	for (count d = 2; d <= dMax; d++) { // FIXME: <= ?
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
			ng_[g] = d + 1;
			if ((b_[g] * (d + 1)) > (ndRest_[d] + nBulk_[d])) {
				// special handling of last group
				if (b_[g] != 1) {
					throw std::runtime_error("?"); // TODO: what happens here? check for indexing error
				}
				ng_[g] = (ndRest_[d] + nBulk_[d]);
			}
			DEBUG("rhoStar: " << rhoStar);
			assert (d <= dMax);
			DEBUG("c_[d] : " << c_[d]);
			rhoStar = std::pow(c_[d], (1.0 / 3.0));
			dStar = (ng_[g] - 1) * rhoStar;
			assert (dStar <= d);
			wBulk_[d] = 0.5 * nBulk_[d] * (d - dStar);
			assert (wBulk_[d] > 0);
			DEBUG("rhoStar: " << rhoStar);
			assert ( (1 - rhoStar) > 0); // avoid division by 0
			wg_[g] = b_[g] * 0.5 * ng_[g] * (ng_[g] - 1) * std::log(1.0 / (1 - rhoStar)); // correct log?
			nFillStar = (b_[g] * ng_[g]) - nBulk_[d];
		} else {
			wBulk_[d] = 0;
		}
		wd_[d] = wFill_[d] + wBulk_[d];
		r_[d] = wFill_[d] / wd_[d];
	}

	// DEBUG: print all arrays for debugging

	DEBUG("nd_ degree distribution: nd_[d] = number of nodes with degree d : " << Aux::vectorToString(nd_));
	DEBUG("c_ clustering coefficient per degree: " << Aux::vectorToString(c_));
	DEBUG("id_ : index i_d for first node of each degree d" << Aux::vectorToString(id_));
	DEBUG("nFill_  number of filler nodes per degree " << Aux::vectorToString(nFill_));
	DEBUG("wd_ sum of wFill and wBulk: " << Aux::vectorToString(wd_));
	DEBUG("r_ ratio of fill excess degree for degree : " << Aux::vectorToString(r_));
	DEBUG("ig_ start index for affinity group : " << Aux::vectorToString(ig_));
	DEBUG("b_ b_[g]: number of blocks in a group : " << Aux::vectorToString(b_));
	DEBUG("wg_ wg_[g]: weight of the group g : " << Aux::vectorToString(wg_));
	DEBUG("ng_  ng_[g]: number of blocks in the affinity group g : " << Aux::vectorToString(ng_));

	DEBUG("nBulk_ number of bulk nodes of degree d : " << Aux::vectorToString(nBulk_));
	DEBUG("ndRest_ : " << Aux::vectorToString(ndRest_));
	DEBUG("wFill_ weight of the degree-d fill nodes for phase 2 : " << Aux::vectorToString(wFill_));
	DEBUG("wBulk_  : weight of degree-d bulk nodes for phase 2 " << Aux::vectorToString(wBulk_));

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


	DEBUG("start sampling");
	// TODO: parallelize this
	for (count j = 0; j < w; ++j) {
		double r = this->rand.probability();
		if (r < (w1 / w)) {
			this->samplePhaseOne();
		} else {
			this->samplePhaseTwo();
		}
	}

}

void BTERGenerator::samplePhaseOne() {
	index g;
	do {
		g = this->rand.choice(wd_); // choose group
	} while (g == none);
	// FIXME: when g = 0 u and v become garbage
	double r1 = this->rand.probability();
	double r2 = this->rand.probability();
	double r3 = this->rand.probability();
	index delta = ig_[g] + std::floor(r1 * b_[g]) * ng_[g]; // choose block and compute its offset  // FIXME: invalid read
	assert (delta >= 0);
	node u = std::floor(r2 * ng_[g]) + delta;
	node v = std::floor(r3 * (ng_[g] - 1)) + delta;
	if (v == u) { // TODO: why? should this be (v >= u) as in the pseudocode?
		v += 1;
	}
	TRACE("adding phase 1 edge (" << u << ", " << v << ")");
	this->G->addEdge(u, v); // TODO: decouple edge insertion from edge generation for parallelism - alternative: make graph thread-safe
}

void BTERGenerator::samplePhaseTwo() {
	node u = this->samplePhaseTwoNode();
	node v = this->samplePhaseTwoNode();
	TRACE("adding phase 2 edge (" << u << ", " << v << ")");
	this->G->addEdge(u, v); // TODO: decouple edge insertion from edge generation for parallelism
}

std::vector<count> BTERGenerator::generatePowerLawDegreeDistribution(count n, double gamma, double a) {

	// power-law function
	auto f = [&](double x){
		return a * std::pow(x, -1.0 * gamma);
	};

	// sample the function and create sum for normalization
	double s = 0.0;
	for (index i = 1; i < n; ++i) {
		s += f(i);
	}

	std::vector<count> dist(n, 0); // the degree distribution
	// sample again and calculate node counts per degree
	for (index i = 1; i < n; ++i) {
		dist[i] = (f(i) / s) * n;
	}

	return dist;
}

Clustering BTERGenerator::getAffinityBlocks() {

	// TODO: get block start and end indices

}

node BTERGenerator::samplePhaseTwoNode() {
	degree d = none;
	do {
		d = this->rand.choice(wd_);
	} while (d == none); // first entry in wd_ is none (-1)
	double r1 = this->rand.probability();
	double r2 = this->rand.probability();
	node u;
	count nFill = nFill_[d];
	assert (nFill >= 0);
	if (r1 < nFill) {
		u = std::floor(r2 * nFill) + id_[d];
	} else {
		u = std::floor(r2 * (nd_[d] - nFill)) + (id_[d] + nFill);
	}

	assert (u >= 0);
	return u;
}

} /* namespace NetworKit */
