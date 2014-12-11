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

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator() : nodes(0), currentfactor(0), alpha(1), stretch(1), moveEachStep(0), factorgrowth(0), moveDistance(0) {
	initializeQuadTree();
	initializeMovement();
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
	initializeQuadTree();
	initializeMovement();
	R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(nodes);
	r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
}

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(vector<double> &angles, vector<double> &radii, double stretch, double initialFactor, double moveEachStep, double factorgrowth, double moveDistance) {
	this->angles = angles;
	this->radii = radii;
	this->nodes = angles.size();
	this->stretch = stretch;
	assert(radii.size() == nodes);
	R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(nodes);
	r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	quad = Quadtree<index>(r);
	currentfactor = initialFactor;
	this->alpha = 1;//not needed any more
	this->moveEachStep = moveEachStep;
	this->factorgrowth = factorgrowth;
	this->moveDistance = moveDistance;
	for (index i = 0; i < nodes; i++) {
		assert(radii[i] < r);
		quad.addContent(i, angles[i], radii[i]);
	}
	this->initialized = true;
	INFO("Filled Quadtree");
	initializeMovement();
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
	double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(nodes);
	angles.resize(nodes);
	radii.resize(nodes);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	quad = Quadtree<index>(r);
	HyperbolicSpace::fillPoints(&angles, &radii, stretch, alpha);
	INFO("Generated Points");
	for (index i = 0; i < nodes; i++) {
		assert(radii[i] < R);
		quad.addContent(i, angles[i], radii[i]);
	}
	INFO("Filled Quadtree");
}

Graph DynamicHyperbolicGenerator::getGraph() const {
	double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(nodes);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	/**
	 * The next call is unnecessarily expensive, since it constructs a new QuadTree.
	 * Reduces code duplication, though.
	 */
	return HyperbolicGenerator().generate(angles, radii, r, currentfactor*R);
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

	for (index step = 0; step < nSteps; step++) {
		if (factorgrowth != 0) {
			getEventsFromFactorGrowth(result);
		}

		if (moveEachStep > 0 && moveDistance > 0) {
			getEventsFromNodeMovement(result);
		}
		result.push_back(GraphEvent(GraphEvent::TIME_STEP));
	}
	return result;
}

void DynamicHyperbolicGenerator::getEventsFromFactorGrowth(vector<GraphEvent> &result) {
	if (currentfactor == 0 && factorgrowth < 0) {
		return;
	}
	double newfactor = currentfactor + factorgrowth;
	if (newfactor < 0) newfactor = 0;

	double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(nodes);

	/**
	 * TODO: get all neighbours in the beginning, sort them by hyperbolic distance, move along edge array.
	 */
	#pragma omp parallel for
	for (index i = 0; i < nodes; i++) {
		assert(R*newfactor > R*currentfactor);
		vector<index> oldset = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*currentfactor);
		//we only add new edges, don't remove any. The order of the points should be the same
		vector<index> newset = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*newfactor);
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

void DynamicHyperbolicGenerator::moveNode(index toMove) {
	double hyperbolicRadius = HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[toMove]);

	//angular movement
	double maxcdf = cosh(alpha*R);
	double mincdf = 1;
	double currcdf = cosh(alpha*hyperbolicRadius);

	double newcosh = currcdf + alpha*radialMovement[toMove];
	double newphi = angles[toMove];
	//bounce off the boundary
	if (newcosh > maxcdf) {
		newcosh -= 2*(newcosh - maxcdf);
		radialMovement[toMove] *= -1;
		DEBUG("Node ", toMove, " bounced off upper boundary, radial movement set to ", radialMovement[toMove]);
	}
	if (newcosh < mincdf) {
		newcosh += 2*(mincdf - newcosh);
		radialMovement[toMove] *= -1;
		DEBUG("Node ", toMove, " crossed center, radial movement set to ", radialMovement[toMove]);

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
	newphi += angularMovement[toMove]/newradius;
	if (newphi < 0) newphi += (floor(-newphi/(2*M_PI))+1)*2*M_PI;
	if (newphi > 2*M_PI) newphi -= floor(newphi/(2*M_PI))*2*M_PI;

	newradius = HyperbolicSpace::hyperbolicRadiusToEuclidean(newradius);
	if (newradius >= r) newradius = std::nextafter(newradius, std::numeric_limits<double>::lowest());

	angles[toMove] = newphi;
	radii[toMove] = newradius;
}

void DynamicHyperbolicGenerator::getEventsFromNodeMovement(vector<GraphEvent> &result) {
	count oldStreamMarker = result.size();
	vector<index> toWiggle;
	vector<vector<index> > oldNeighbours;
	//TODO: One could parallelize this.
	for (index i = 0; i < nodes; i++) {
		if (Aux::Random::real(1) < moveEachStep) {
			toWiggle.push_back(i);
			oldNeighbours.push_back(quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), R*currentfactor));
		}
	}
	/**
	 * Tried to parallelize this, didn't bring any benefit.
	 * Not surprising, since most of the work - manipulating the QuadTree - needs to be done in a critical section
	 */
	#pragma omp parallel for
	for (index j = 0; j < toWiggle.size(); j++) {
		//wiggle this node!

		double oldphi = angles[toWiggle[j]];
		double oldr = radii[toWiggle[j]];
		moveNode(toWiggle[j]);
		//updating Quadtree
		#pragma omp critical
		{
			bool removed = quad.removeContent(toWiggle[j], oldphi, oldr);
			assert(removed);
			quad.addContent(toWiggle[j], angles[toWiggle[j]], radii[toWiggle[j]]);
		}
	}

	//now get the new edges and see what changed
	#pragma omp parallel for
	for (index j = 0; j < toWiggle.size(); j++) {
		vector<index> newNeighbours = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[toWiggle[j]], radii[toWiggle[j]]), R*currentfactor);
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

} /* namespace NetworKit */
