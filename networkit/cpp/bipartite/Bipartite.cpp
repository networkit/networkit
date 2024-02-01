/*
 * Bipartite.cpp
 *
 * Created on: 18.09.2023
 *     Author: Michael Kaibel
 */

#include <networkit/bipartite/Bipartite.hpp>
#include <deque>

namespace NetworKit {

Bipartite::Bipartite(const Graph &G) : G(&G) {
    partition = Partition(G.upperNodeIdBound(), none);
}

bool Bipartite::isBipartite() {
    assureFinished();

    return bipartite;
}

const Partition &Bipartite::getPartition() {
    assureFinished();

    if (!bipartite)
        throw std::runtime_error("Can't provide bipartition on non-bipartite graph");

    return partition;
}

const std::vector<node> &Bipartite::getOddCycle() {
    assureFinished();

    if (bipartite)
        throw std::runtime_error("Can't provide odd circle on bipartit graph");

    return oddCircle;
}

void Bipartite::run() {
    hasRun = true;

    std::vector<node> parent(G->upperNodeIdBound(), none);

    for (node v : G->nodeRange()) {
        if (partition[v] != none)
            continue;

        partition[v] = 0;

        std::deque queue(1, v);

        while (!queue.empty()) {
            node w = queue.front();
            queue.pop_front();

            assert(partition[w] != none);

            for (node x : G->neighborRange(w)) {
                if (partition[w] == partition[x]) {
                    bipartite = false;
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

    bipartite = true;
}

void Bipartite::findOddCircle(std::vector<node> &parent, NetworKit::node v, NetworKit::node w) {
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

} // namespace NetworKit
