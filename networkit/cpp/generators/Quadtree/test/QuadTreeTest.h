/*
 * QuadTreeTest.h
 *
 *  Created on: 28.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef QUADTREETEST_H_
#define QUADTREETEST_H_

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "../Quadtree.h"

using std::vector;

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity

class QuadTreeTest: public testing::Test {
public:
	QuadTreeTest();
	virtual ~QuadTreeTest();

protected:
	template <class T>
	QuadNode<T> getRoot(Quadtree<T> &tree) {
		return tree.root;
	}

	template <class T>
	vector<QuadNode<T> > getChildren(QuadNode<T> &node) {
		return node.children;
	}
};

} /* namespace NetworKit */
#endif /* QUADTREETEST_H_ */
