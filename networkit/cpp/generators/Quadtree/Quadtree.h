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
	Quadtree(double maxR,bool theoreticalSplit=false, double alpha=1, count capacity=1000, bool diagnostics = false) {
		root = QuadNode<T>(0, 0, 2*M_PI, maxR, capacity, 0,theoreticalSplit,alpha,diagnostics);
		this->maxRadius = maxR;
	}

	Quadtree(count n, double stretch, bool theoreticalSplit=false, double alpha=1, count capacity=1000, bool diagnostics = false) {
		double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
		count numberOfThreads = omp_get_num_threads();
		count k = ceil(log(numberOfThreads)/log(4));
		root = QuadNode<T>(0, 0, 2*M_PI, R, capacity, 0,theoreticalSplit,alpha,diagnostics);
		fillParallel(n, alpha, k, 0, root);
		//bernoulli to distribute
	}


	void fillInParallel(count l, double alpha, count seqDepth, count offset, QuadNode<T> &currentNode) {
		if (seqDepth > 0) {
			if (currentNode.isLeaf) currentNode.split();
			for (int i = 0; i < currentNode.children.size(); i++) {
				fillInParallel(l/4, alpha, seqDepth -1, offset+(l/4)*i, currentNode.children[i]);
			}
		} else {
				#pragma omp task
				{
					vector<double> angles(l);
					vector<double> radii(l);
					HyperbolicSpace::fillPoints(angles, radii, currentNode.leftAngle, currentNode.rightAngle, currentNode.minR, currentNode.maxR, alpha);
					for (index i = 0; i < l; i++) {
						currentNode.addContent(i+offset, angles[i], radii[i]);
					}
				}
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

	/**
	 * Get elements whose hyperbolic distance to the query point is less than the hyperbolic distance
	 *
	 * Safe to call in parallel if diagnostics were not activated
	 *
	 * @param circleCenter Cartesian coordinates of the query circle's center
	 * @param hyperbolicRadius Radius of the query circle
	 */
	vector<T> getElementsInHyperbolicCircle(Point2D<double> circleCenter, double hyperbolicRadius) {
		Point2D<double> origin(0,0);
		vector<T> circleDenizens;
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

		/**
		 * return the elements
		 */
		return circleDenizens;
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

	/**
	 * trims the vectors used to hold the content in the leaf cells. Reduces memory usage, makes changes slower
	 */
	void trim() {
		root.trim();
	}

	/**
	 * Reset the counters of necessary and unnecessary comparisons.
	 */
	void resetCounter() {
		root.resetCounter();
	}

private:
	QuadNode<T> root;
	double maxRadius;
};
}

#endif /* QUADTREE_H_ */
