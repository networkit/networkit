/*
 * GenericDijkstra.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#ifndef GENERICDIJKSTRA_H_
#define GENERICDIJKSTRA_H_

#include "Graph.h"
#include "APP.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Generic Dijkstra single source algorithm.
 */
template <class T>
class GenericDijkstra : public APP<T> {

// friend class DynDijkstra;
// friend class DynDijkstra2;

public:

	/**
	 * Creates the Dijkstra class for @a G and the source node @a source.
	 *
	 * @param G The graph.
	 * @param source The source node.
	 * @param storePaths	store paths and number of paths?
	 */
    GenericDijkstra(const Graph& G, node source, std::vector<T> edgeWeights, bool storePaths=true, bool storeStack=false, node target = none) : APP<T>(G, source, storePaths, storeStack, target, edgeWeights) {

    }

	/**
	 * Performs the Dijkstra SSSP algorithm on the graph given in the constructor.
	 */
    void run() {

        TRACE("initializing Dijkstra data structures");
        // init distances
        T infDist = T::getZero();

        distances.clear();
        distances.resize(G.upperNodeIdBound(), infDist);
        if (storePaths) {
            previous.clear();
            previous.resize(G.upperNodeIdBound());
            npaths.clear();
            npaths.resize(G.upperNodeIdBound(), 0);
            npaths[source] = 1;
        }

        if (storeStack) {
            std::vector<node> empty;
            std::swap(stack, empty);
        }
        distances[source] = T::getOne();
        // priority queue with distance-node pairs
        distances[source] = T::getOne();

        Aux::PrioQueue<T, node> pq(distances);

        auto relax([&](node u, node v, edgeid e) {
            T dist = distances[v] + T(distances[u] * edgeWeights[e]);
            if (dist != distances[v]) {
                distances[v] = T(dist);
                if (storePaths) {
                    previous[v] = {u}; // new predecessor on shortest path
                    npaths[v] = npaths[u];
                }
                TRACE("Decreasing key of ", v);
                TRACE("pq size: ", pq.size());
                pq.decreaseKey(distances[v], v);
                TRACE("pq size: ", pq.size());
            } else if (storePaths) {
                previous[v].push_back(u); 	// additional predecessor
                npaths[v] += npaths[u]; 	// all the shortest paths to u are also shortest paths to v now
            }
        });

        bool breakWhenFound = (target != none);
        TRACE("traversing graph");
        while (pq.size() > 0) {
            TRACE("pq size: ", pq.size());
            node current = pq.extractMin().second;
            TRACE("current node in Dijkstra: " , current);
            TRACE("pq size: ", pq.size());
            if (breakWhenFound && target == current) {
                break;
            }

            if (storeStack) {
                stack.push_back(current);
            }

            G.forEdgesOf(current, relax);
        }

    }

private:

    using APP<T>::G;
    using APP<T>::source;
    using APP<T>::target;
    using APP<T>::distances;
    using APP<T>::previous;
    using APP<T>::npaths;
    using APP<T>::stack;
    using APP<T>::storePaths;
    using APP<T>::storeStack;
    using APP<T>::edgeWeights;
};

} /* namespace NetworKit */
#endif /* GENERICDIJKSTRA_H_ */
