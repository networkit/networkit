/*
 * BFSTree.h
 *
 *  Created on: Nov 21, 2013
 *      Author: gbrueckner
 */

#ifndef BFSTREE_H_
#define BFSTREE_H_

#include "../graph/Graph.h"
#include <unordered_map>

namespace NetworKit {

/**
 * TODO: class documentation
 */
class BFSTree : public Graph {
    std::unordered_map<node, count> bfs;
    node _deepest;
    count _depth;
    bool _spanning;
public:
	BFSTree(const Graph& G, node root);
	virtual ~BFSTree();
    count depth();

    // Returns a node with maximum depth in the tree.
    node deepest();

    // Returns true iff this tree is a spanning tree.
    bool spanning();
};

} /* namespace NetworKit */
#endif /* BFSTREE_H_ */
