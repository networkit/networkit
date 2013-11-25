/*
 * BFSTree.cpp
 *
 *  Created on: Nov 21, 2013
 *      Author: gbrueckner 
 */

#include "BFSTree.h"
#include <queue>
#include <utility>

namespace NetworKit {

BFSTree::BFSTree(const Graph& G, node root) : Graph(G.numberOfNodes()) {

    using namespace std;

    count infDist = numeric_limits<count>::max();

    _depth = 0;
    _deepest = root;
    _spanning = true;

    G.forNodes([&](node v) {
            bfs[v] = infDist;
        });
    bfs[root] = 0;

    std::queue<node> q;
    q.push(root);

    while (!q.empty()) {
        node u = q.front();
        q.pop();
        G.forNeighborsOfInRandomOrder(u, [&](node v) {
                if (bfs[v] == infDist) {
                    bfs[v] = bfs[u] + 1;
                    _depth = max(_depth, bfs[v]);
                    if (bfs[v] == _depth)
                        _deepest = v;
                    q.push(v);
                    this->addEdge(u, v);
                }
            });
    }

    this->forNodes([&](node v) {
        if (degree(v) < 1)
            this->removeNode(v);
    });

    G.forNodes([&]() {
            return _spanning;
        },
        [&](node v) {
            if (bfs[v] == infDist)
                _spanning = false;
        });
}

BFSTree::~BFSTree() {

}

bool BFSTree::spanning() {
    return _spanning;
}

count BFSTree::depth() {
    return _depth;
}

node BFSTree::deepest() {
    return _deepest;
}

} /* namespace NetworKit */
