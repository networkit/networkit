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
#include "QuadNode.h"
#include "../HyperbolicSpace.h"

namespace NetworKit {

template <class T>
class Quadtree {
	friend class QuadTreeTest;
public:
	Quadtree() {
		root = QuadNode<T>();
	}

	Quadtree(double maxR) {
		root = QuadNode<T>(0, 0, 2*M_PI, maxR, 20, 0);
	}

	virtual ~Quadtree() {

	}

	void addContent(T newcomer, double angle, double R) {
		root.addContent(newcomer, angle, R);
	}
	vector<T> getElements() {
		return root.getElements();
	}

	vector<T> getCloseElements(Point<double> query, double maxDistance) {
		Point<double> origin(0,0);
		vector<T> circleDenizens;
		vector<Point<double> > positions;
		double hyperbolicFromOrigin = HyperbolicSpace::getHyperbolicDistance(origin, query);
		if (hyperbolicFromOrigin < maxDistance) {
			/*
			 * circle will overlap origin, approach not feasible!
			 */
			double maxR = pow((cosh(maxDistance+hyperbolicFromOrigin)-1)/(cosh(maxDistance+hyperbolicFromOrigin)+1), 0.5);
			vector<T> subresult;
			vector<Point<double> > subpos;
			root.getElementsInEuclideanCircle(0, 2*M_PI, 0, maxR, origin, maxR, subresult, subpos);
			assert(subresult.size() == subpos.size());
			//filter manually. Sigh.
			DEBUG("Filter manually with radius ", maxR);
			for (index i = 0; i < subresult.size(); i++) {
				if (HyperbolicSpace::getHyperbolicDistance(subpos[i], query) < maxDistance) {
					circleDenizens.push_back(subresult[i]);
				}
			}
		} else {
		/**
		 * compute circle and polar rectangle
		 */
		Point<double> pointOnEdge = HyperbolicSpace::getPointOnHyperbolicCircle(query, maxDistance);
		double distance = HyperbolicSpace::getHyperbolicDistance(query, pointOnEdge);
		assert(abs(distance - maxDistance) < 0.00001);
		Point<double> center;
		double radius, minPhi, maxPhi;
		HyperbolicSpace::getEuclideanCircle(query, pointOnEdge, center, radius);
		DEBUG("Using circle at (", center[0], ",",center[1], ") with radius ", radius);
		double minR = center.length() - radius;
		double maxR = center.length() + radius;
		assert(maxR < 1);
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
			 * what to do if they overlap the 2\pi line? Well, have to make two separate calls and collect
			 */
		}

		/**
		 * get Elements in Circle
		 */


		root.getElementsInEuclideanCircle(minPhi, maxPhi, minR, maxR, center, radius, circleDenizens, positions);
		if (minPhi < 0) {
			root.getElementsInEuclideanCircle(2*M_PI+minPhi, 2*M_PI, minR, maxR, center, radius, circleDenizens, positions);
		}
		if (maxPhi > 2*M_PI) {
			root.getElementsInEuclideanCircle(0, maxPhi - 2*M_PI, minR, maxR, center, radius, circleDenizens, positions);
		}

		}

		/**
		 * return them
		 */
		return circleDenizens;
		//return root.getCloseElements(query, maxDistance);
	}


private:
	QuadNode<T> root;
};
}

#endif /* QUADTREE_H_ */
