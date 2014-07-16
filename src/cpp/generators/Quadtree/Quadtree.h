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
		Point<double> center;
		double hyperbolicFromOrigin = HyperbolicSpace::getHyperbolicDistance(origin, query);
		Point<double> pointOnEdge;
		if (hyperbolicFromOrigin < maxDistance) {
			/*
			 * TODO: move code to method
			 */
			//use right triangle at origin to construct new point at phi_p, r_p
			double phi_q, r_q;
			HyperbolicSpace::cartesianToPolar(query, phi_q, r_q);
			double phi_p, r_p;
			phi_p = phi_q + M_PI/2;
			if (phi_p > 2*M_PI) phi_p -= 2*M_PI;
			double r_p_nom = (cosh(maxDistance)-1)*(1-r_q*r_q)-2*r_q*r_q;
			double r_p_denom = (2+(cosh(maxDistance)-1)*(1-r_q*r_q));
			r_p = pow(r_p_nom/r_p_denom, 0.5);
			pointOnEdge = HyperbolicSpace::polarToCartesian(phi_p, r_p);
		} else {
			/**
		 	 * compute circle and polar rectangle
		 	 */
			pointOnEdge = HyperbolicSpace::getPointOnHyperbolicCircle(query, maxDistance);
		}
		double distance = HyperbolicSpace::getHyperbolicDistance(query, pointOnEdge);
		assert(abs(distance - maxDistance) < 0.00001);

		double minPhi, maxPhi, radius;
		HyperbolicSpace::getEuclideanCircle(query, pointOnEdge, center, radius);
		TRACE("Using circle at (", center[0], ",",center[1], ") with radius ", radius);
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


		root.getElementsInEuclideanCircle(minPhi, maxPhi, minR, maxR, center, radius, circleDenizens);
		if (minPhi < 0) {
			root.getElementsInEuclideanCircle(2*M_PI+minPhi, 2*M_PI, minR, maxR, center, radius, circleDenizens);
		}
		if (maxPhi > 2*M_PI) {
			root.getElementsInEuclideanCircle(0, maxPhi - 2*M_PI, minR, maxR, center, radius, circleDenizens);
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
