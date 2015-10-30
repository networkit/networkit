/*
 * RandomSpanningTree.h
 *
 *  Created on: 20.06.2015
 *      Author: Henning
 */

#ifndef RANDOMSPANNINGTREE_H_
#define RANDOMSPANNINGTREE_H_

#include "../graph/Graph.h"

namespace NetworKit {

class RandomSpanningTree {
public:
	RandomSpanningTree(const Graph& G);
	virtual ~RandomSpanningTree();

	void run();

	Graph getTree();

private:
	const Graph& g;
	Graph tree;
};

} /* namespace NetworKit */
#endif /* RANDOMSPANNINGTREE_H_ */
