/*
 * DynamicHyperbolicGenerator.cpp
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#include <cmath>

#include "DynamicHyperbolicGenerator.h"
#include "RHGGenerator.h"
#include "../geometric/HyperbolicSpace.h"
#include "../auxiliary/Parallel.h"

using std::vector;
namespace NetworKit {

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(count n, double avgDegree, double exp, double moveEachStep, double moveDistance) {
	nodeCount = n;
	this->alpha = (exp-1)/2;
	this->moveEachStep = moveEachStep;
	this->moveDistance = moveDistance;
	this->initialized = false;
	R = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha);
	initializePoints();
	initializeMovement();
	initializeBands();
}

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(std::vector<double> &angles, std::vector<double> &radii, double R, double alpha, double moveEachStep, double moveDistance) {
	this->angles = angles;
	this->radii = radii;
	this->nodeCount = angles.size();
	this->alpha = alpha;
	assert(radii.size() == nodeCount);
	this->R = R;
	for (double r : radii) {
		assert(r < R);
	}

	this->moveEachStep = moveEachStep;
	this->moveDistance = moveDistance;

	this->initialized = true;

	initializeMovement();
	initializeBands();
}

void DynamicHyperbolicGenerator::initializeMovement() {
	/**
	 *
	 */
	angularMovement.resize(nodeCount);
	radialMovement.resize(nodeCount);
	int scale = 10;
	for (index i = 0; i < nodeCount; i++) {
		angularMovement[i] = Aux::Random::real(-moveDistance, moveDistance);
		radialMovement[i] = Aux::Random::real(-scale*moveDistance, scale*moveDistance);
	}
	INFO("Initialized Movement");
}

void DynamicHyperbolicGenerator::initializePoints() {
	assert(nodeCount > 0);
	if (initialized) return;
	else initialized = true;
	angles.resize(nodeCount);
	radii.resize(nodeCount);
	RHGGenerator::fillPoints(angles, radii, R, alpha);

	INFO("Generated Points");
}

void DynamicHyperbolicGenerator::initializeBands() {
	//1.Generate bandRadius'
	bandRadii = RHGGenerator::getBandRadii(nodeCount, R);
	INFO("Got Band Radii");
	assert(bandRadii.size() > 1);
	//2. Initialize empty bands
	bands.clear();
	bands.resize(bandRadii.size() - 1);
	bandAngles.resize(bandRadii.size() - 1);
	assert(angles.size() == nodeCount);
	assert(radii.size() == nodeCount);
	//3. Put points to bands
	INFO("Starting Point distribution");
	#pragma omp parallel for
	for (index j = 0; j < bands.size(); j++){
		for (index i = 0; i < nodeCount; i++){
			if (radii[i] >= bandRadii[j] && radii[i] <= bandRadii[j+1]){
				bands[j].push_back(Point2D<double>(angles[i], radii[i], i));
				bandAngles[j].push_back(angles[i]);
			}
		}
	}
	INFO("Filled Bands");
}

Graph DynamicHyperbolicGenerator::getGraph() const {
	/**
	 * The next call is unnecessarily expensive, since it constructs a new QuadTree.
	 * Reduces code duplication, though.
	 */
	return RHGGenerator().generate(angles, radii, R);
}

std::vector<Point<float> > DynamicHyperbolicGenerator::getCoordinates() const {
	const count n = angles.size();
	assert(radii.size() == n);
	std::vector<Point<float> > result;
	for (index i = 0; i < n; i++) {
		Point2D<double> coord = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
		Point<float> temp(coord[0], coord[1]);
		result.push_back(temp);
	}
	return result;
}

std::vector<GraphEvent> DynamicHyperbolicGenerator::generate(count nSteps) {
	if (!initialized) {
		assert(angles.size() == 0);
		assert(radii.size() == 0);
		initializePoints();
		initializeMovement();
		initializeBands();
	}
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
	double hyperbolicRadius = radii[toMove];

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

	angles[toMove] = newphi;
	radii[toMove] = newradius;
}

vector<index> DynamicHyperbolicGenerator::getNeighborsInBands(index i, bool bothDirections) {
	const double phi = angles[i];
	const double r = radii[i];
	assert(bands.size() == bandAngles.size());
	assert(bands.size() == bandRadii.size() -1);
	count expectedDegree = (4/M_PI)*nodeCount*exp(-(radii[i])/2);
	vector<index> near;
	near.reserve(expectedDegree*1.1);
	for(index j = 0; j < bands.size(); j++){
		if(bothDirections || bandRadii[j+1] > radii[i]){
			double minTheta, maxTheta;
			std::tie (minTheta, maxTheta) = RHGGenerator::getMinMaxTheta(phi, r, bandRadii[j], R);

			vector<Point2D<double>> neighborCandidates = RHGGenerator::getPointsWithinAngles(minTheta, maxTheta, bands[j], bandAngles[j]);

			const count sSize = neighborCandidates.size();
			for(index w = 0; w < sSize; w++){
				if(HyperbolicSpace::nativeHyperbolicDistance(phi, r, neighborCandidates[w].getX(), neighborCandidates[w].getY()) <= R) {
					const index possibleNeighbor = neighborCandidates[w].getIndex();
					if(possibleNeighbor != i){
						near.push_back(possibleNeighbor);
					}
				}
			}
		}
	}
	return near;
}

void DynamicHyperbolicGenerator::getEventsFromNodeMovement(vector<GraphEvent> &result) {
	for (index b = 0; b < bandAngles.size(); b++) {
		assert(std::is_sorted(bandAngles[b].begin(), bandAngles[b].end()));
	}

	count oldStreamMarker = result.size();
	vector<index> toWiggle;
	vector<vector<index> > oldNeighbours;
	//TODO: One could parallelize this.
	for (index i = 0; i < nodeCount; i++) {
		if (Aux::Random::real(1) < moveEachStep) {

			toWiggle.push_back(i);
			Point2D<double> q = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
			vector<index> localOldNeighbors = getNeighborsInBands(i, true);

			oldNeighbours.push_back(localOldNeighbors);
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

		int bandIndex = -1;
		for (index j = 0; j < bands.size(); j++){
			if ((oldr) < bandRadii[j+1]) {
				bandIndex = oldr;
				break;
			}
		}
		assert(bandIndex >= 0);



		//TODO: find old and new band of wiggled node, find old and new positions
		//moveNode(toWiggle[j]);
		//updating Quadtree
		#pragma omp critical
		{
			//bool removed = quad.removeContent(toWiggle[j], oldphi, oldr);
			//assert(removed);
			//quad.addContent(toWiggle[j], angles[toWiggle[j]], radii[toWiggle[j]]);
		}
	}

	//now get the new edges and see what changed
	#pragma omp parallel for
	for (index j = 0; j < toWiggle.size(); j++) {
		vector<index> newNeighbours = getNeighborsInBands(toWiggle[j], true);;

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
	Aux::Parallel::sort(result.begin()+oldStreamMarker, result.end(), GraphEvent::compare);
	auto end = std::unique(result.begin()+oldStreamMarker, result.end(), GraphEvent::equal);
	result.erase(end, result.end());

}

} /* namespace NetworKit */
