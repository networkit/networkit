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
#include "../auxiliary/Parallel.h"

using std::vector;
namespace NetworKit {

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(count n, double avgDegree, double exp, double moveEachStep, double moveDistance) {
	nodes = n;
	this->alpha = (exp-1)/2;
	this->moveEachStep = moveEachStep;
	this->moveDistance = moveDistance;
	this->initialized = false;
	R = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha, 0);
	r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	initializeQuadTree();
	initializeMovement();
}
DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(std::vector<double> &angles, std::vector<double> &radii, double avgDegree, double exp, double moveEachStep, double moveDistance) {
	this->angles = angles;
	this->radii = radii;
	this->nodes = angles.size();
	this->alpha = (exp-1)/2;
	assert(radii.size() == nodes);
	R = HyperbolicSpace::getTargetRadius(nodes, nodes*avgDegree/2, alpha, 0);
	r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	quad = Quadtree<index>(r);
	this->moveEachStep = moveEachStep;
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
	angles.resize(nodes);
	radii.resize(nodes);
	quad = Quadtree<index>(r);
	double oldR = HyperbolicSpace::hyperbolicAreaToRadius(nodes);
	HyperbolicSpace::fillPoints(angles, radii, R / oldR, alpha);
	INFO("Generated Points");
	for (index i = 0; i < nodes; i++) {
		assert(radii[i] < R);
		quad.addContent(i, angles[i], radii[i]);
	}
	INFO("Filled Quadtree");
}

Graph DynamicHyperbolicGenerator::getGraph() const {
	/**
	 * The next call is unnecessarily expensive, since it constructs a new QuadTree.
	 * Reduces code duplication, though.
	 */
	return HyperbolicGenerator().generate(angles, radii, r, R);
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

		if (moveEachStep > 0 && moveDistance > 0) {
			getEventsFromNodeMovement(result);
		}
		result.push_back(GraphEvent(GraphEvent::TIME_STEP));
	}
	return result;
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
	bool suppressLeft = false;

	count oldStreamMarker = result.size();
	vector<index> toWiggle;
	vector<vector<index> > oldNeighbours;
	//TODO: One could parallelize this.
	for (index i = 0; i < nodes; i++) {
		if (Aux::Random::real(1) < moveEachStep) {
			vector<index> localOldNeighbors;
			toWiggle.push_back(i);
			Point2D<double> q = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
			quad.getElementsInHyperbolicCircle(q, R, suppressLeft, localOldNeighbors);
			oldNeighbours.push_back(localOldNeighbors);//add temperature
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
		vector<index> newNeighbours;
		Point2D<double> q = HyperbolicSpace::polarToCartesian(angles[toWiggle[j]], radii[toWiggle[j]]);
		quad.getElementsInHyperbolicCircle(q, R, suppressLeft, newNeighbours);

		Aux::Parallel::sort(oldNeighbours[j].begin(), oldNeighbours[j].end());
		Aux::Parallel::sort(newNeighbours.begin(), newNeighbours.end());
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
	Aux::Parallel::sort(result.begin()+oldStreamMarker, result.end(), GraphEvent::compare);
	auto end = std::unique(result.begin()+oldStreamMarker, result.end(), GraphEvent::equal);
	result.erase(end, result.end());
}

} /* namespace NetworKit */
