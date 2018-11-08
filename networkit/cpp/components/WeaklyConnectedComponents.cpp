/*
 * WeaklyConnectedComponents.cpp
 *
 *  Created on: June 20, 2017
 *      Author: Eugenio Angriman
*/

#include "WeaklyConnectedComponents.h"

namespace NetworKit {

    WeaklyConnectedComponents::WeaklyConnectedComponents(const Graph& G) :
    G(G) {
        if (!G.isDirected()) {
            throw std::runtime_error("Weakly Connected Components can be computeed for directed graphs. Use ConnectedComponents for undirected graphs.");
        }
    }

    void WeaklyConnectedComponents::init() {

        if (hasRun) {
            compSize.clear();
            std::fill(components.begin(), components.end(), none);
        }
        else {
            components.assign(G.upperNodeIdBound(), none);
        }

        hasRun = false;
    }

    void WeaklyConnectedComponents::run() {

        // Initialization of data structures
        init();

        // Queue for BFS.
        std::queue<node> q;

        // Perform BFSs to assign a component ID to each node.
        G.forNodes([&](node u) {

            // Node u has not been visited.
            if (components[u] == none) {

                // New component ID.
                index c = compSize.size();
                components[u] = c;
                compSize.insert(std::make_pair(c, 1));


                // Start a new BFS from u.
                q.push(u);

                do {
                    node v = q.front();
                    q.pop();

                    // Enqueue neighbors (both from in and out edges) and set
                    // new component.
                    G.forNeighborsOf(v, [&](node w) {
                        updateComponent(c, w, q, false);
                    });

                    G.forInNeighborsOf(v, [&](node w) {
                        updateComponent(c, w, q, true);

                    });
                } while (!q.empty());
            }
        });

        hasRun = true;
    }

    void WeaklyConnectedComponents::updateComponent(
        index c, node w, std::queue<node>& q, bool inNeighbor
    ) {

        if (components[w] == none) {
            q.push(w);
            components[w] = c;
            compSize.find(c)->second += 1;
        }
    }


    std::vector<std::vector<node> > WeaklyConnectedComponents::getComponents() {
        assureFinished();

        // transform partition into vector of unordered_set
        std::vector<std::vector<node> > result(compSize.size());
        std::map<index, count> compIndex;

        int i = 0;
        for (auto it=compSize.begin(); it!=compSize.end(); ++it) {
            auto indexIterator = compIndex.find(it->first);
            if (indexIterator == compIndex.end()) {
                compIndex.insert(std::pair<index, count>(it->first, i));
                ++i;
            }
        }

        G.forNodes([&](node u) {
            result[compIndex.find(components[u])->second].push_back(u);
        });

        return result;
    }
}
