/*
 * Quadtree.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <vector>
#include <memory>
#include <cmath>
#include <omp.h>
#include <functional>
#include "QuadNode.h"
#include "../../geometric/HyperbolicSpace.h"
#include "../../auxiliary/Parallel.h"

namespace NetworKit {

template <class T, bool poincare=true>
class Quadtree {
	friend class QuadTreeGTest;
public:
	Quadtree() {
		root = QuadNode<T, poincare>();
		this->maxRadius = 1;
	}

	/**
	 * @param maxR Radius of the managed area. Must be smaller than 1.
	 * @param theoreticalSplit If true, split cells to get the same area in each child cell. Default is false
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param capacity How many points can inhabit a leaf cell before it is split up?
	 *
	 */
	Quadtree(double maxR,bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		root = QuadNode<T,poincare>(0, 0, 2*M_PI, maxR, capacity, theoreticalSplit,alpha,balance);
		this->maxRadius = maxR;
	}

	count fillInParallel(count l, double alpha, count seqThreshold, count offset, QuadNode<T, poincare> &currentNode) {
		if (l > seqThreshold) {
			if (currentNode.height() == 1) currentNode.split();
			double treeArea = HyperbolicSpace::effectiveAreaInCell(currentNode.getLeftAngle(), currentNode.getRightAngle(), currentNode.getMinR(), currentNode.getMaxR(), alpha);
			count coveredPoints = 0;
			double coveredArea = 0;
			for (index i = 0; i < currentNode.children.size(); i++) {
				count pointsInSubtree;
				if (i < currentNode.children.size()-1) {
					double subTreeArea = HyperbolicSpace::effectiveAreaInCell(currentNode.children[i].getLeftAngle(), currentNode.children[i].getRightAngle(), currentNode.children[i].getMinR(), currentNode.children[i].getMaxR(), alpha);
					assert(treeArea-coveredArea >= subTreeArea);
					std::binomial_distribution<int> distribution(l-coveredPoints,subTreeArea/(treeArea-coveredArea));
					pointsInSubtree = distribution(Aux::Random::getURNG());
					coveredArea += subTreeArea;
				} else {
					pointsInSubtree = l-coveredPoints;
				}
				coveredPoints += pointsInSubtree;
				offset = fillInParallel(pointsInSubtree, alpha, seqThreshold, offset, currentNode.children[i]);
				//offset += pointsInSubtree;
			}
		} else {
				#pragma omp task shared(currentNode) firstprivate(offset)
				{
					vector<double> angles(l);
					vector<double> radii(l);
					HyperbolicSpace::fillPoints(angles, radii, currentNode.getLeftAngle(), currentNode.getRightAngle(),							HyperbolicSpace::EuclideanRadiusToHyperbolic(currentNode.getMinR()),
							HyperbolicSpace::EuclideanRadiusToHyperbolic(currentNode.getMaxR()), alpha);
					for (index i = 0; i < l; i++) {
						currentNode.addContent(i+offset, angles[i], HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]));
					}
				}
				offset += l;
			}
			return offset;
		}

	Quadtree(count n, double stretch, bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		static_assert(poincare, "Parallel Quadtree filling only works with Poincaré disk");
		double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
		double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
		count numberOfThreads = omp_get_max_threads();
		//double k = ceil(log(numberOfThreads)/log(4));
		root = QuadNode<T,poincare>(0, 0, 2*M_PI, r, capacity, theoreticalSplit,alpha,balance);
		maxRadius = r;
		count result;
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				result = fillInParallel(n, alpha, n/numberOfThreads, 0, root);
			}
		}
		assert(result == n);
		root.recount();
		assert(root.size() == n);
	}

	Quadtree(const vector<double> &angles, const vector<double> &radii, const vector<T> &content, double stretch, bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		const count n = angles.size();
		assert(angles.size() == radii.size());
		assert(radii.size() == content.size());
		double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
		double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
		root = QuadNode<T>(0, 0, 2*M_PI, r, capacity, theoreticalSplit,alpha,balance);
		maxRadius = r;
		for (index i = 0; i < n; i++) {
			assert(content[i] < n);
			root.addContent(content[i], angles[i], radii[i]);
		}
	}

	/**
	 * @param newcomer content to be added at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	void addContent(T newcomer, double angle, double r) {
		root.addContent(newcomer, angle, r);
	}

	/**
	 * @param newcomer content to be removed at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	bool removeContent(T toRemove, double angle, double r) {
		return root.removeContent(toRemove, angle, r);
	}

	/**
	 * Get all elements, regardless of position
	 *
	 * @return vector<T> of elements
	 */
	vector<T> getElements() const {
		return root.getElements();
	}

	void extractCoordinates(vector<double> &anglesContainer, vector<double> &radiiContainer) const {
		root.getCoordinates(anglesContainer, radiiContainer);
	}

	/**
	 * Get elements whose hyperbolic distance to the query point is less than the hyperbolic distance
	 *
	 *
	 * @param circleCenter Cartesian coordinates of the query circle's center
	 * @param hyperbolicRadius Radius of the query circle
	 */
	vector<T> getElementsInHyperbolicCircle(Point2D<double> circleCenter, double hyperbolicRadius) const {
		vector<T> circleDenizens;
		getElementsInHyperbolicCircle(circleCenter, hyperbolicRadius, circleDenizens);
		return circleDenizens;
	}

	void getElementsInHyperbolicCircle(const Point2D<double> circleCenter, const double hyperbolicRadius, const bool suppressLeft, vector<T> &circleDenizens) const {
		assert(circleDenizens.empty());
		double cc_phi, cc_r;
		HyperbolicSpace::cartesianToPolar(circleCenter, cc_phi, cc_r);
		//Transform hyperbolic circle into Euclidean circle
		double minPhi, maxPhi, radius, r_e;
		HyperbolicSpace::getEuclideanCircle(cc_r, hyperbolicRadius, r_e, radius);
		Point2D<double> center = HyperbolicSpace::polarToCartesian(cc_phi, r_e);
		double minR = r_e - radius;
		double maxR = r_e + radius;
		//assert(maxR < 1);//this looks fishy
		if (maxR > 1) maxR = 1;
		if (minR < 0) {
			maxR = std::max(abs(minR), maxR);
			minR = 0;
			minPhi = 0;
			maxPhi = 2*M_PI;
		} else {
			double spread = asin(radius / r_e);
			//double phi_c, r_c;
			//HyperbolicSpace::cartesianToPolar(center, phi_c, r_c);
			minPhi = cc_phi - spread;
			maxPhi = cc_phi + spread;
			/**
			 * If the circle overlaps the 2\pi line, we have to make two separate calls and collect
			 */
		}

		if (suppressLeft) minPhi = cc_phi;
		/**
		 * get Elements in Euclidean circle
		 */

		bool wraparound = false;
		root.getElementsInEuclideanCircle(center, radius, circleDenizens, minPhi, maxPhi, minR, maxR);
		if (minPhi < 0) {
			root.getElementsInEuclideanCircle(center, radius, circleDenizens, 2*M_PI+minPhi, 2*M_PI, minR, maxR);
			wraparound = true;
		}
		if (maxPhi > 2*M_PI) {
			root.getElementsInEuclideanCircle(center, radius, circleDenizens, 0, maxPhi - 2*M_PI, minR, maxR);
			wraparound = true;
		}

		for (T denizen : circleDenizens) {
			if (denizen >= size()) {
				DEBUG("Content ", denizen, " found in quadtree of size ", size(), ", as one of ", circleDenizens.size(), " neighbours.");
			}
			assert(denizen < size());//TODO: remove this after debugging, in general the quadtree should handle arbitrary contents
		}

		//we have sort(deg(v)) here! This is not good, but does not make the asymptotical complexity of O(deg(v) log n) worse.
		if (wraparound) {
			Aux::Parallel::sort(circleDenizens.begin(), circleDenizens.end());
			auto newend = unique(circleDenizens.begin(), circleDenizens.end());
			count toRemove = circleDenizens.end() - newend;
			count remaining = newend - circleDenizens.begin();
			if (toRemove > 0) {
				DEBUG("Removing, ", toRemove, " duplicate entries, keeping ", remaining);
				circleDenizens.resize(remaining);
			}
		}

		for (T denizen : circleDenizens) {
			if (denizen >= size()) DEBUG("Content ", denizen, " found in quadtree of size ", size(), ", as one of ", circleDenizens.size(), " neighbours, after sorting");
			assert(denizen < size());//TODO: remove this after debugging, in general the quadtree should handle arbitrary contents
		}
	}

	void getElementsInHyperbolicCircle(const Point2D<double> circleCenter, const double hyperbolicRadius, vector<T> &circleDenizens) const {
		getElementsInHyperbolicCircle(circleCenter, hyperbolicRadius, false, circleDenizens);
	}

	count getElementsProbabilistically(Point2D<double> euQuery, std::function<double(double)> prob, vector<T> &circleDenizens) {
		return root.getElementsProbabilistically(euQuery, prob, false, circleDenizens);
	}

	count getElementsProbabilistically(Point2D<double> euQuery, std::function<double(double)> prob, bool suppressLeft, vector<T> &circleDenizens) {
		return root.getElementsProbabilistically(euQuery, prob, suppressLeft, circleDenizens);
	}

	void recount() {
		root.recount();
	}

	count size() const {
		return root.size();
	}

	count height() const {
		return root.height();
	}

	count countLeaves() const {
		return root.countLeaves();
	}

	index indexSubtree(index nextID) {
		return root.indexSubtree(nextID);
	}

	index getCellID(double phi, double r) const {
		return root.getCellID(phi, r);
	}

	double getMaxRadius() const {
		return maxRadius;
	}

	void sortPointsInLeaves() {
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				root.sortPointsInLeaves();
			}
		}
	}

	void reindex() {
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				root.reindex(0);
			}
		}
	}

	/**
	 * trims the vectors used to hold the content in the leaf cells. Reduces memory usage, makes changes slower
	 */
	void trim() {
		root.trim();
	}

private:
	QuadNode<T, poincare> root;
	double maxRadius;
};
}

#endif /* QUADTREE_H_ */
