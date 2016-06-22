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

DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(count n, double avgDegree, double exp, double T, double moveEachStep, double moveDistance) {
	nodeCount = n;
	this->alpha = (exp-1)/2;
	this->T = T;
	this->moveEachStep = moveEachStep;
	this->moveDistance = moveDistance;
	this->initialized = false;
	R = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha, T);
	initializePoints();
	initializeMovement();
	if (T > 0) {
		initializeQuadTree();
	} else {
		recomputeBands();
	}
}
DynamicHyperbolicGenerator::DynamicHyperbolicGenerator(std::vector<double> &angles, std::vector<double> &radii, double R, double alpha, double T, double moveEachStep, double moveDistance) {
	this->angles = angles;
	this->radii = radii;
	this->nodeCount = angles.size();
	this->alpha = alpha;
	this->T = T;

	assert(radii.size() == nodeCount);
	this->R = R;
	for (double r : radii) {
		assert(r < R);
	}

	this->moveEachStep = moveEachStep;
	this->moveDistance = moveDistance;

	this->initialized = true;

	initializeMovement();
	if (T > 0) {
		initializeQuadTree();
	} else {
		recomputeBands();
	}
}

void DynamicHyperbolicGenerator::initializePoints() {
	assert(nodeCount > 0);
	if (initialized) return;
	else initialized = true;
	angles.resize(nodeCount);
	radii.resize(nodeCount);
	HyperbolicSpace::fillPoints(angles, radii, R, alpha);

	INFO("Generated Points");
}

void DynamicHyperbolicGenerator::initializeMovement() {
	angularMovement.resize(nodeCount);
	radialMovement.resize(nodeCount);
	int scale = 10;
	for (index i = 0; i < nodeCount; i++) {
		angularMovement[i] = Aux::Random::real(-moveDistance, moveDistance);
		radialMovement[i] = Aux::Random::real(-scale*moveDistance, scale*moveDistance);
	}
}

void DynamicHyperbolicGenerator::initializeQuadTree() {
	for (index i = 0; i < nodeCount; i++) {
		assert(radii[i] < R);
		quad.addContent(i, angles[i], radii[i]);
	}
	INFO("Filled Quadtree");
}

void DynamicHyperbolicGenerator::recomputeBands() {
	//1.Generate bandRadius'
	bandRadii = HyperbolicGenerator::getBandRadii(nodeCount, R);
	INFO("Got Band Radii");
	assert(bandRadii.size() > 1);
	//2. Initialize empty bands
	bands.clear();
	bands.resize(bandRadii.size() - 1);
	bandAngles.clear();
	bandAngles.resize(bandRadii.size() - 1);
	assert(angles.size() == nodeCount);
	assert(radii.size() == nodeCount);

	//ensure points are sorted
	vector<index> permutation(nodeCount);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	//3. Put points to bands
	INFO("Starting Point distribution");
	#pragma omp parallel for
	for (index j = 0; j < bands.size(); j++){
		for (index i = 0; i < nodeCount; i++){
			double alias = permutation[i];
			if (radii[alias] >= bandRadii[j] && radii[alias] <= bandRadii[j+1]){
				bands[j].push_back(Point2D<double>(angles[alias], radii[alias], alias));
				bandAngles[j].push_back(angles[alias]);
			}
		}
	}
	INFO("Filled Bands");
}

Graph DynamicHyperbolicGenerator::getGraph() const {
	/**
	 * The next call is unnecessarily expensive, since it constructs a new QuadTree / bands.
	 * Reduces code duplication, though.
	 */
	return HyperbolicGenerator().generate(angles, radii, R, T);
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
		if (T > 0) {
			initializeQuadTree();
		} else {
			recomputeBands();
		}
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
	const double r = radii[i];
	const double phi = angles[i];
	//const double coshr = cosh(radii[i]);
	//const double sinhr = sinh(radii[i]);
	//const double coshR = cosh(R);
	assert(bands.size() == bandAngles.size());
	assert(bands.size() == bandRadii.size() -1);
	count expectedDegree = (4/M_PI)*nodeCount*exp(-(radii[i])/2);
	vector<index> near;
	near.reserve(expectedDegree*1.1);
	for(index j = 0; j < bands.size(); j++){
		if(bothDirections || bandRadii[j+1] > radii[i]){
			double minTheta, maxTheta;
			std::tie (minTheta, maxTheta) = HyperbolicGenerator::getMinMaxTheta(phi, r, bandRadii[j], R);

			vector<Point2D<double>> neighborCandidates = HyperbolicGenerator::getPointsWithinAngles(minTheta, maxTheta, bands[j], bandAngles[j]);

			const count sSize = neighborCandidates.size();
			for(index w = 0; w < sSize; w++){
				if (HyperbolicSpace::nativeDistance(phi, r, neighborCandidates[w].getX(), neighborCandidates[w].getY()) <= R) {
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
	bool suppressLeft = false;

	//now define lambda
	double tresholdDistance = R;
	double beta = 1/T;
	//dummy lambda:
	std::function<double(double)> edgeProb = [beta, tresholdDistance](double distance) -> double {return distance <= tresholdDistance ? 1 : 0;};

	if (T > 0) {
		if (!(beta == beta)) {
			DEBUG("Value of beta is ", beta, ", which is invalid. T=", T);
		}
		assert(T == 0 || beta == beta);

		edgeProb = [beta, tresholdDistance](double distance) -> double {return 1 / (exp(beta*(distance-tresholdDistance)/2)+1);};
	}

	count oldStreamMarker = result.size();
	vector<index> toWiggle;
	vector<vector<index> > oldNeighbours;
	//TODO: One could parallelize this.
	for (index i = 0; i < nodeCount; i++) {
		if (Aux::Random::real(1) < moveEachStep) {
			vector<index> localOldNeighbors;
			toWiggle.push_back(i);
			if (T == 0) {
				localOldNeighbors = getNeighborsInBands(i, true);
			} else {
				Point2D<double> q = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
				quad.getElementsProbabilistically(q, edgeProb, suppressLeft, localOldNeighbors);
			}
			oldNeighbours.push_back(localOldNeighbors);
		}
	}

	#pragma omp parallel for
	for (index j = 0; j < toWiggle.size(); j++) {
		//wiggle this node!

		double oldphi = angles[toWiggle[j]];
		double oldr = radii[toWiggle[j]];
		moveNode(toWiggle[j]);
		if (T > 0) {
			//updating Quadtree
			#pragma omp critical
			{
				bool removed = quad.removeContent(toWiggle[j], oldphi, oldr);
				assert(removed);
				quad.addContent(toWiggle[j], angles[toWiggle[j]], radii[toWiggle[j]]);
			}
		}
	}

	if (T == 0) {
		recomputeBands();//update bands
	}

	//now get the new edges and see what changed
	#pragma omp parallel for
	for (index j = 0; j < toWiggle.size(); j++) {
		vector<index> newNeighbours;
		if (T == 0) {
			newNeighbours = getNeighborsInBands(toWiggle[j], true);
		} else {
			Point2D<double> q = HyperbolicSpace::polarToCartesian(angles[toWiggle[j]], radii[toWiggle[j]]);
			quad.getElementsProbabilistically(q, edgeProb, suppressLeft, newNeighbours);
		}

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
