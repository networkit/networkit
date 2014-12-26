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
#include "QuadNode.h"
#include "../../geometric/HyperbolicSpace.h"

namespace NetworKit {

template <class T>
class Quadtree {
	friend class QuadTreeTest;
public:
	Quadtree() {
		root = QuadNode<T>();
		this->maxRadius = 1;
	}

	/**
	 * @param maxR Radius of the managed area. Must be smaller than 1.
	 * @param theoreticalSplit If true, split cells to get the same area in each child cell. Default is false
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param capacity How many points can inhabit a leaf cell before it is split up?
	 * @param diagnostics Count how many necessary and unnecessary comparisons happen in leaf cells? Will cause race condition and false sharing in parallel calls
	 *
	 */
	Quadtree(double maxR,bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		root = QuadNode<T>(0, 0, 2*M_PI, maxR, capacity, 0,theoreticalSplit,alpha,balance);
		this->maxRadius = maxR;
	}

	count fillInParallel(count l, double alpha, count seqThreshold, count offset, QuadNode<T> &currentNode) {
		if (l > seqThreshold) {
			if (currentNode.height() == 1) currentNode.split();
			double treeArea = HyperbolicSpace::effectiveAreaInCell(currentNode.getLeftAngle(), currentNode.getRightAngle(), currentNode.getMinR(), currentNode.getMaxR(), alpha);
			for (int i = 0; i < currentNode.children.size(); i++) {
				double subTreeArea = HyperbolicSpace::effectiveAreaInCell(currentNode.children[i].getLeftAngle(), currentNode.children[i].getRightAngle(), currentNode.children[i].getMinR(), currentNode.children[i].getMaxR(), alpha);
				const count pointsInSubtree = l*(subTreeArea/treeArea);
				offset = fillInParallel(pointsInSubtree, alpha, seqThreshold, offset, currentNode.children[i]);
				//offset += pointsInSubtree;
			}
		} else {
				#pragma omp task shared(currentNode) firstprivate(offset)
				{
					vector<double> angles(l);
					vector<double> radii(l);
					HyperbolicSpace::fillPoints(angles, radii, currentNode.getLeftAngle(), currentNode.getRightAngle(),
							HyperbolicSpace::EuclideanRadiusToHyperbolic(currentNode.getMinR()),
							HyperbolicSpace::EuclideanRadiusToHyperbolic(currentNode.getMaxR()), alpha);
					for (index i = 0; i < l; i++) {
						currentNode.addContent(i+offset, angles[i], radii[i]);
					}
				}
				offset += l;
			}
			return offset;
		}

	Quadtree(count n, double stretch, bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
		double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
		count numberOfThreads = omp_get_max_threads();
		//double k = ceil(log(numberOfThreads)/log(4));
		root = QuadNode<T>(0, 0, 2*M_PI, r, capacity, 0,theoreticalSplit,alpha,balance);
		count result;
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				result = fillInParallel(n, alpha, n/numberOfThreads, 0, root);
			}
		}
		vector<double> missingAngles(n-result);
		vector<double> missingRadii(n-result);
		HyperbolicSpace::fillPoints(missingAngles, missingRadii, stretch, alpha);
		for (index i = result; i < n; i++) {
			root.addContent(i, missingAngles[i-result], missingRadii[i-result]);
		}
		//bernoulli to distribute
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
	 * Safe to call in parallel if diagnostics were not activated
	 *
	 * @param circleCenter Cartesian coordinates of the query circle's center
	 * @param hyperbolicRadius Radius of the query circle
	 */
	vector<T> getElementsInHyperbolicCircle(Point2D<double> circleCenter, double hyperbolicRadius) {
		vector<T> circleDenizens;
		getElementsInHyperbolicCircle(circleCenter, hyperbolicRadius, circleDenizens);
		return circleDenizens;
	}

	void getElementsInHyperbolicCircle(Point2D<double> circleCenter, double hyperbolicRadius, vector<T> &circleDenizens) {
		Point2D<double> origin(0,0);
		Point2D<double> center;

		//Transform hyperbolic circle into Euclidean circle
		double minPhi, maxPhi, radius;
		HyperbolicSpace::getEuclideanCircle(circleCenter, hyperbolicRadius, center, radius);
		double minR = center.length() - radius;
		double maxR = center.length() + radius;
		//assert(maxR < 1);//this looks fishy
		if (maxR > 1) maxR = 1;
		if (minR < 0) {
			maxR = std::max(abs(minR), maxR);
			minR = 0;
			minPhi = 0;
			maxPhi = 2*M_PI;
		} else {
			double spread = asin(radius / center.length());
			double phi_c, r_c;
			HyperbolicSpace::cartesianToPolar(center, phi_c, r_c);
			minPhi = phi_c - spread;
			maxPhi = phi_c + spread;
			/**
			 * If the circle overlaps the 2\pi line, we have to make two separate calls and collect
			 */
		}

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

		//we have sort(deg(v)) here! This is not good, but does not make the asymptotical complexity of O(deg(v) log n) worse.
		if (wraparound) {
			std::sort(circleDenizens.begin(), circleDenizens.end());
			auto newend = unique(circleDenizens.begin(), circleDenizens.end());
			circleDenizens.resize(newend - circleDenizens.begin());
		}
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

	index getCellID(double phi, double r) {
		return root.getCellID(phi, r);
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
	QuadNode<T> root;
	double maxRadius;
};
}

#endif /* QUADTREE_H_ */
