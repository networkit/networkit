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

namespace NetworKit {

template <class T>
class Quadtree {
public:
	Quadtree() {

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
	vector<T> getCloseElements(double angle, double R, double maxDistance) {
		return root.getCloseElements(angle, R, maxDistance);
	}

private:
	QuadNode<T> root;
};
}

#endif /* QUADTREE_H_ */
