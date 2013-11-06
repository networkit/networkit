/*
 * ConductanceTree.h
 *
 *  Created on: Aug 9, 2013
 *      Author: Henning
 */

#ifndef CONDUCTANCETREE_H_
#define CONDUCTANCETREE_H_

#include "../graph/Graph.h"
#include "../clustering/Clustering.h"
#include <set>

namespace NetworKit {

class ConductanceTree {
public:
	ConductanceTree();
	virtual ~ConductanceTree();

	Clustering bestCutInBfsTree(const Graph& g, node root);
};

} /* namespace NetworKit */
#endif /* CONDUCTANCETREE_H_ */
