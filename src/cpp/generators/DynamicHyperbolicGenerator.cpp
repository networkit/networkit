/*
 * DynamicHyperbolicGenerator.cpp
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#include "DynamicHyperbolicGenerator.h"

#include "../geometric/HyperbolicSpace.h"

using std::vector;
namespace NetworKit {

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator() {
	// TODO Auto-generated constructor stub

}

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(count n, double initialFactor, double alpha, double stretch, double moveEachStep, double factorgrowth, double moveDistance) {
	nodes = n;
	currentfactor = initialFactor;
	this->alpha = alpha;
	this->stretch = stretch;
	this->moveEachStep = moveEachStep;
	this->factorgrowth = factorgrowth;
	this->moveDistance = moveDistance;
	this->initialized = false;
}

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(vector<double> &angles, vector<double> &radii, double R, double initialFactor, double moveEachStep, double factorgrowth, double moveDistance) {


	this->angles = angles;
	this->radii = radii;
	this->nodes = angles.size();
	assert(radii.size() == nodes);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	quad = Quadtree<index>(r);
	currentfactor = initialFactor;
	this->alpha = 1;//not needed any more
	this->stretch = 1;//not needed any more
	this->moveEachStep = moveEachStep;
	this->factorgrowth = factorgrowth;
	this->moveDistance = moveDistance;
	this->initialized = true;
	for (index i = 0; i < nodes; i++) {
		assert(radii[i] < r);
		quad.addContent(i, angles[i], radii[i]);
	}
	INFO("Filled Quadtree");
}

DynamicHyperbolicGenerator::~DynamicHyperbolicGenerator() {
	// TODO Auto-generated destructor stub
}

void DynamicHyperbolicGenerator::initializeGraph() {
	if (initialized) return;
	else initialized = true;
	double R = stretch*acosh((double)nodes/(2*M_PI)+1);
	angles.resize(nodes);
	radii.resize(nodes);
	double rad_nom = (cosh(R)-1);
	double rad_denom = (cosh(R)+1);
	double r = sqrt(rad_nom/rad_denom);
	quad = Quadtree<index>(r);
	HyperbolicSpace::fillPoints(&angles, &radii, stretch, alpha);
	INFO("Generated Points");
	for (index i = 0; i < nodes; i++) {
		assert(radii[i] < R);
		quad.addContent(i, angles[i], radii[i]);
	}
	INFO("Filled Quadtree");
}

std::vector<GraphEvent> DynamicHyperbolicGenerator::generate(count nSteps) {
	if (!initialized) initializeGraph();
	assert(quad.size() == nodes);
	vector<GraphEvent> result;
	double R = stretch*acosh((double)nodes/(2*M_PI)+1);

	for (index step = 0; step < nSteps; step++) {
		assert(factorgrowth == 0 || moveEachStep == 0 || moveDistance == 0);
		if (factorgrowth != 0) {
			//nodes are stationary, growing factors
			double newfactor = currentfactor + factorgrowth;
/**
 * most efficient way: get all neighbours in the beginning, sort them by hyperbolic distance, move along edge array
 */
			for (index i = 0; i < nodes; i++) {
				vector<index> oldset = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*currentfactor);
				//we only add new edges, don't remove any. The order of the points should be the same
				vector<index> newset = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*newfactor);
				assert(newset.size() >= oldset.size());

				//these should not be necessary, the sets should already be in the same order. This is suspicious.
				std::sort(oldset.begin(), oldset.end());
				std::sort(newset.begin(), newset.end());
				index oldindex = 0;
				index newindex = 0;
				for (newindex = 0; newindex < newset.size(); newindex++) {
					double distance = HyperbolicSpace::getHyperbolicDistance(angles[i], radii[i], angles[newset[newindex]], radii[newset[newindex]]);
					if (oldindex < oldset.size() && newset[newindex] == oldset[oldindex]) {
						//skip element
						assert(distance <= R*currentfactor);
						TRACE("Skipping old edge (", i, ", ", newset[newindex], ")");
						oldindex++;
					} else if (i < newset[newindex]){
						assert(distance <= R*newfactor);
						assert(distance >= R*currentfactor);
						result.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, i, newset[newindex]));
					}
				}
				assert(oldindex == oldset.size());
			}
			currentfactor = newfactor;
		}
		else {
			vector<index> toWiggle;
			vector<vector<index> > oldNeighbours;
			for (index i = 0; i < nodes; i++) {
				if (Aux::Random::real(1) < moveEachStep) {
					toWiggle.push_back(i);
					oldNeighbours.push_back(quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*currentfactor));
				}
			}
			for (index j = 0; j < toWiggle.size(); j++) {
				//wiggle this node!
				double hyperbolicRadius = HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[toWiggle[j]]);
				Point2D<double> point = HyperbolicSpace::polarToCartesian(angles[toWiggle[j]], hyperbolicRadius);
				Point2D<double> offset = HyperbolicSpace::polarToCartesian(Aux::Random::real(2*M_PI), moveEachStep);
				double newphi, newradius;
				HyperbolicSpace::cartesianToPolar(point + offset, newphi, newradius);

				bool removed = quad.removeContent(toWiggle[j], angles[toWiggle[j]], radii[toWiggle[j]]);
				assert(removed);

				angles[toWiggle[j]] = newphi;
				radii[toWiggle[j]] = newradius;
				quad.addContent(toWiggle[j], newphi, newradius);
			}

			//now get the new edges and see what changed
			for (index j = 0; j < toWiggle.size(); j++) {
				vector<index> newNeighbours = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[toWiggle[j]], radii[toWiggle[j]]), R*currentfactor);
				std::sort(oldNeighbours[j].begin(), oldNeighbours[j].end());
				std::sort(newNeighbours.begin(), newNeighbours.end());
				index aindex = 0;
				index bindex = 0;
				while (aindex < oldNeighbours[j].size() && bindex < newNeighbours.size()) {
					if (oldNeighbours[j][aindex] == newNeighbours[bindex]) {
						//unchanged edge, skip
						aindex++;
						bindex++;
					} else if (oldNeighbours[j][aindex] < newNeighbours[bindex]) {
						//edge only in one set.
						result.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, toWiggle[j], oldNeighbours[j][aindex]));
						aindex++;
					} else {
						result.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, toWiggle[j], newNeighbours[bindex]));
						bindex++;
					}
				}

				//handling leftover edges

				while (aindex < oldNeighbours[j].size()) {
					result.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, toWiggle[j], oldNeighbours[j][aindex]));
					aindex++;
				}

				while (bindex < newNeighbours.size()) {
					result.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, toWiggle[j], newNeighbours[bindex]));
					bindex++;
				}
			}

		}
		result.push_back(GraphEvent(GraphEvent::TIME_STEP));
	}
	return result;
}



} /* namespace NetworKit */
