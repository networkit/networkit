/*
* Bipartit.cpp
*
* Created on: 18.09.2023
*     Author: Michael Kaibel
*/

#include <networkit/bipartit/Bipartit.hpp>
#include <deque>

namespace NetworKit {

Bipartit::Bipartit(const Graph &G) : G(&G) {
    partition = Partition(G.upperNodeIdBound(), none);
}

bool Bipartit::isBipartit() {
    assureFinished();

    return bipartit;
}

const Partition &Bipartit::getPartition() {
    assureFinished();

    if (not bipartit)
        throw std::runtime_error("Can't provide bipartition on non-bipartite graph");

    return partition;
}

const std::vector<node> &Bipartit::getOddCircle() {
    assureFinished();

    if (bipartit)
        throw std::runtime_error("Can't provide odd circle on bipartit graph");

    return oddCircle;
}

void Bipartit::run() {
    hasRun = true;

    std::vector<node> parent(G->upperNodeIdBound(), none);

    for (node v : G->nodeRange()) {
        if (partition[v] != none)
            continue;

        partition[v] = 0;

        std::deque queue(1, v);

        while (not queue.empty()) {
            node w = queue.front();
            queue.pop_front();

            assert(partition[w] != none);

            for (node x : G->neighborRange(w)) {
                if (partition[w] == partition[x]) {
                    bipartit = false;
                    findOddCircle(parent, w, x);
                    return;
                }

                if (partition[x] == none) {
                    queue.emplace_back(x);
                    partition[x] = 1 - partition[w];
                    parent[x] = w;
                }
            }
        }
    }

    bipartit = true;
}

void Bipartit::findOddCircle(std::vector<node> &parent, NetworKit::node v, NetworKit::node w) {
    std::vector<node> pathToW;
    while (v != w) {
        oddCircle.emplace_back(v);
        pathToW.emplace_back(w);
        v = parent[v];
        w = parent[w];
    }

    oddCircle.emplace_back(v);

    for (count i = 0; i < pathToW.size(); i++) {
        oddCircle.emplace_back(pathToW[pathToW.size() - i - 1]);
    }
}

} // NetworKit