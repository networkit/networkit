/*
 * DynamicHyperbolicGenerator.cpp
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#include <cmath>

#include "DynamicHyperbolicGenerator.h"
#include "HyperbolicGenerator.h"
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
	initializeMovement();
}

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(vector<double> &angles, vector<double> &radii, double stretch, double initialFactor, double moveEachStep, double factorgrowth, double moveDistance) {
	this->angles = angles;
	this->radii = radii;
	this->nodes = angles.size();
	this->stretch = stretch;
	assert(radii.size() == nodes);
	double R = stretch*acosh((double)nodes/(2*M_PI)+1);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	quad = Quadtree<index>(r);
	currentfactor = initialFactor;
	this->alpha = 1;//not needed any more
	this->moveEachStep = moveEachStep;
	this->factorgrowth = factorgrowth;
	this->moveDistance = moveDistance;
	this->initialized = true;
	for (index i = 0; i < nodes; i++) {
		assert(radii[i] < r);
		quad.addContent(i, angles[i], radii[i]);
	}
	INFO("Filled Quadtree");
	initializeQuadTree();
	initializeMovement();
}

DynamicHyperbolicGenerator::~DynamicHyperbolicGenerator() {
	// TODO Auto-generated destructor stub
}

void DynamicHyperbolicGenerator::initializeMovement() {
	angularMovement.resize(nodes);
	radialMovement.resize(nodes);
	int scale = 10;
	for (index i = 0; i < nodes; i++) {
		angularMovement[i] = Aux::Random::real(-moveDistance, moveDistance);
		radialMovement[i] = Aux::Random::real(-scale*moveDistance, scale*moveDistance);
	}
}

void DynamicHyperbolicGenerator::initializeQuadTree() {
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

Graph DynamicHyperbolicGenerator::getGraph() {
	/**
	 * while the Quadtree doesn't strictly need to be initialized this early, it has to be done anyway and there is no harm doing it now
	 */
	if (!initialized) initializeQuadTree();
	double R = stretch*acosh((double)nodes/(2*M_PI)+1);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	/**
	 * The next call is unnecessarily expensive, since it constructs a new QuadTree.
	 * Reduces code duplication, though.
	 */
	return HyperbolicGenerator::generate(&angles, &radii, r, currentfactor*R);
}

std::vector<Point<float> > DynamicHyperbolicGenerator::getCoordinates() const {
	count n = angles.size();
	assert(radii.size() == n);
	std::vector<Point<float> > result;
	for (index i = 0; i < angles.size(); i++) {
		Point2D<double> coord = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
		Point<float> temp(coord[0], coord[1]);
		result.push_back(temp);
	}
	return result;
}

std::vector<Point<float> > DynamicHyperbolicGenerator::getHyperbolicCoordinates() const {
	count n = angles.size();
	assert(radii.size() == n);
	std::vector<Point<float> > result;
	for (index i = 0; i < angles.size(); i++) {
		Point2D<double> coord = HyperbolicSpace::polarToCartesian(angles[i], HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[i]));
		Point<float> temp(coord[0], coord[1]);
		result.push_back(temp);
	}
	return result;
}

std::vector<GraphEvent> DynamicHyperbolicGenerator::generate(count nSteps) {
	if (!initialized) initializeQuadTree();
	assert(quad.size() == nodes);
	vector<GraphEvent> result;
	double R = stretch*acosh((double)nodes/(2*M_PI)+1);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);

	for (index step = 0; step < nSteps; step++) {
		count oldStreamMarker = result.size();
		assert(factorgrowth == 0 || moveEachStep == 0 || moveDistance == 0);
		if (factorgrowth != 0) {
			//nodes are stationary, growing factors
			double newfactor = currentfactor + factorgrowth;
/**
 * TODO: get all neighbours in the beginning, sort them by hyperbolic distance, move along edge array.
 */
			#pragma omp parallel for
			for (index i = 0; i < nodes; i++) {
				assert(R*newfactor > R*currentfactor);
				vector<index> oldset = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*currentfactor);
				//we only add new edges, don't remove any. The order of the points should be the same
				vector<index> newset = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*newfactor);
				assert(newset.size() >= oldset.size());
				std::sort(oldset.begin(), oldset.end());
				std::sort(newset.begin(), newset.end());
				vector<index> difference(newset.size());

				//get new edges
				auto it = std::set_difference(newset.begin(), newset.end(), oldset.begin(), oldset.end(), difference.begin());

				//keep only those pointing to higher node indices, we don't want any multiedges
				it = std::remove_if(difference.begin(), it, [i](index edge){return i >= edge;});
				difference.resize(it - difference.begin());

				#pragma omp critical
				{
					for (auto edge : difference) {
						result.emplace_back(GraphEvent::EDGE_ADDITION, i, edge);
					}
				}

			}
			currentfactor = newfactor;
		}
		else {
			vector<index> toWiggle;
			vector<vector<index> > oldNeighbours;
			//TODO: One could parallelize this.
			for (index i = 0; i < nodes; i++) {
				if (Aux::Random::real(1) < moveEachStep) {
					toWiggle.push_back(i);
					oldNeighbours.push_back(quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*currentfactor));
				}
			}
			/**
			 * Tried to parallelize this, didn't bring any benefit.
			 * Not surprising, since most of the work - manipulating the QuadTree - needs to be done in a critical section
			 */
			for (index j = 0; j < toWiggle.size(); j++) {
				//wiggle this node!
				double hyperbolicRadius = HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[toWiggle[j]]);

				//angular movement

				double maxcdf = cosh(R);
				double mincdf = 1;
				double currcdf = cosh(hyperbolicRadius);

				double newcosh = currcdf + radialMovement[toWiggle[j]];
				double newphi = angles[toWiggle[j]];
				//bounce off the boundary
				if (newcosh > maxcdf) {
					newcosh -= 2*(newcosh - maxcdf);
					radialMovement[toWiggle[j]] *= -1;
					DEBUG("Node ", toWiggle[j], " bounced off upper boundary, radial movement set to ", radialMovement[toWiggle[j]]);
				}
				if (newcosh < mincdf) {
					newcosh += 2*(mincdf - newcosh);
					radialMovement[toWiggle[j]] *= -1;
					DEBUG("Node ", toWiggle[j], " bounced off lower boundary, radial movement set to ", radialMovement[toWiggle[j]]);

					/**
					 * nodes should cross the center, not bounce off it
					 */
					if (newphi > M_PI) {
						newphi -= M_PI;
					}
					else {
						newphi += M_PI;
					}
				}
				double newradius = acosh(newcosh)/alpha;
				//assert(abs(newradius - hyperbolicRadius) < moveEachStep);
				if (newradius >= R) newradius = std::nextafter(R, std::numeric_limits<double>::lowest());
				assert(newradius < R);
				assert(newradius >= 0);

				//double angleMovement = Aux::Random::real(-moveDistance/hyperbolicRadius, moveDistance/hyperbolicRadius);
				newphi += angularMovement[toWiggle[j]]/newradius;
				if (newphi < 0) newphi += (floor(-newphi/(2*M_PI))+1)*2*M_PI;
				if (newphi > 2*M_PI) newphi -= floor(newphi/(2*M_PI))*2*M_PI;

				newradius = HyperbolicSpace::hyperbolicRadiusToEuclidean(newradius);
				if (newradius >= r) newradius = std::nextafter(newradius, std::numeric_limits<double>::lowest());

				//updating Quadtree
				bool removed = quad.removeContent(toWiggle[j], angles[toWiggle[j]], radii[toWiggle[j]]);
				assert(removed);
				angles[toWiggle[j]] = newphi;
				radii[toWiggle[j]] = newradius;
				quad.addContent(toWiggle[j], newphi, newradius);

			}

			//now get the new edges and see what changed
			#pragma omp parallel for
			for (index j = 0; j < toWiggle.size(); j++) {
				vector<index> newNeighbours = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles[toWiggle[j]], radii[toWiggle[j]]), R*currentfactor);
				std::sort(oldNeighbours[j].begin(), oldNeighbours[j].end());
				std::sort(newNeighbours.begin(), newNeighbours.end());
				vector<index> newEdges(newNeighbours.size());
				auto it = std::set_difference(newNeighbours.begin(), newNeighbours.end(), oldNeighbours[j].begin(), oldNeighbours[j].end(), newEdges.begin());
				newEdges.erase(it, newEdges.end());//trim empty space

				vector<index> brokenEdges(oldNeighbours[j].size() - (newNeighbours.size() - newEdges.size()));//this should be the number of broken edges
				it = std::set_difference(oldNeighbours[j].begin(), oldNeighbours[j].end(), newNeighbours.begin(), newNeighbours.end(), brokenEdges.begin());
				assert(it == brokenEdges.end());

				#pragma omp critical
				{
					for (index edge : newEdges) {
						result.emplace_back(GraphEvent::EDGE_ADDITION, toWiggle[j], edge);
					}
					for (index edge : brokenEdges) {
						result.emplace_back(GraphEvent::EDGE_REMOVAL, toWiggle[j], edge);
					}
				}
			}

			/**
			 * since we didn't know whether the other endpoint of an edge was wiggled or not, we may
			 * have added some of them twice.
			 * Remove the doubles:
			 */
			for (auto it = result.begin()+oldStreamMarker; it < result.end(); it++) {
				if (it->u > it->v) std::swap(it->u, it->v);
			}
			std::sort(result.begin()+oldStreamMarker, result.end(), GraphEvent::compare);
			auto end = std::unique(result.begin()+oldStreamMarker, result.end(), GraphEvent::equal);
			result.erase(end, result.end());
		}
		result.push_back(GraphEvent(GraphEvent::TIME_STEP));
	}
	return result;
}
} /* namespace NetworKit */
