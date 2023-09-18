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
    return bipartit;
}

const Partition &Bipartit::getPartition() {
    return partition;
}

void Bipartit::run() {
    hasRun = true;

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
                    return;
                }

                if (partition[x] != none) {
                    queue.emplace_back(x);
                    partition[x] = 1 - partition[w];
                }
            }
        }
    }

    bipartit = true;
}

} // NetworKit