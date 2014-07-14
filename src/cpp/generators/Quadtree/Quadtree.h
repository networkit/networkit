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
		/**
		 * compute circle and polar rectangle
		 */
		Point<double> pointOnEdge = HyperbolicSpace::getPointOnHyperbolicCircle(query, maxDistance);
		double distance = HyperbolicSpace::getHyperbolicDistance(query, pointOnEdge);
		assert(distance == maxDistance);
		Point<double> center;
		double radius, minPhi, maxPhi;
		HyperbolicSpace::getEuclideanCircle(query, pointOnEdge, center, radius);
		double minR = center.length() - radius;
		double maxR = center.length() + radius;
		assert(maxR < 1);
		if (minR < 0) {
			maxR = std::max(abs(minR), maxR);
			minR = 0;
			minPhi = 0;
			maxPhi = 2*M_PI;
		} else {
			double spread = atan(radius / center.length());
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

		vector<T> circleDenizens = root.getElementsInEuclideanCircle(minPhi, maxPhi, minR, maxR, center, radius);
		if (minPhi < 0) {
			vector<T> additional = root.getElementsInEuclideanCircle(2*M_PI+minPhi, 2*M_PI, minR, maxR, center, radius);
			circleDenizens.insert(circleDenizens.end(), additional.begin(), additional.end());
		}
		if (maxPhi > 2*M_PI) {
			vector<T> additional = root.getElementsInEuclideanCircle(0, maxPhi - 2*M_PI, minR, maxR, center, radius);
			circleDenizens.insert(circleDenizens.end(), additional.begin(), additional.end());
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
