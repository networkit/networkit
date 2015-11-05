/*
 * PseudoRandomSpanningTree.h
 *
 *  Created on: 20.06.2015
 *      Author: Henning
 */

#ifndef PSEUDORANDOMSPANNINGTREE_H_
#define PSEUDORANDOMSPANNINGTREE_H_

#include "../graph/Graph.h"

namespace NetworKit {

struct MyEdge {
	node from;
	node to;
	edgeweight weight;
	bool operator<(const MyEdge& other) const {
		return this->weight >= other.weight; // FIXME
	}
};

class PseudoRandomSpanningTree {
public:
	PseudoRandomSpanningTree(const Graph& G);
	virtual ~PseudoRandomSpanningTree() = default;

	void runShuffle();

	void run();

	Graph getTree();


private:
	const Graph& g;
	Graph tree;
};

} /* namespace NetworKit */
#endif /* PSEUDORANDOMSPANNINGTREE_H_ */
