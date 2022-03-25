/*
 * AStarGeneral.hpp
 *
 *  Created on: 09.08.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_A_STAR_GENERAL_HPP_
#define NETWORKIT_DISTANCE_A_STAR_GENERAL_HPP_

#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/distance/STSP.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * AStar path-finding algoritm with general Heuristic function.
 */
template <typename Heuristic>
class AStarGeneral final : public STSP {
public:
    /**
     * Creates the AStarGeneral class for the graph @a G, the source node @a
     * source, and the target node @a target using a template parameter as
     * heuristic function.
     *
     * @param G The graph.
     * @param heu Heuristic function that takes a node as input and returns a
     * lower bound of the distance of that node to the target node.
     * @param source The source node.
     * @param target The target node.
     * @param storePred If true, the algorithm will also store the predecessors
     * and reconstruct a shortest path from @a source and @a target.
     */
    AStarGeneral(const Graph &G, Heuristic heu, node source, node target, bool storePred = true)
        : STSP(G, source, target, storePred), heu(heu), heap{prio} {}

    /*
     * Executes the AStar algorithm.
     */
    void run() override {
        if (source == target) {
            distance = 0.;
            hasRun = true;
            return;
        }

        init();
        constexpr auto infdist = std::numeric_limits<edgeweight>::max();
        const count n = G->upperNodeIdBound();
        distance = infdist;

        std::fill(distFromSource.begin(), distFromSource.end(), infdist);
        distFromSource.resize(n, infdist);
        distFromSource[source] = 0.;

        std::fill(prio.begin(), prio.end(), infdist);
        prio.resize(n, infdist);
        prio[source] = 0.;

        heap.clear();
        heap.reserve(n);
        heap.push(source);

        node top = none;
        do {
            top = heap.extract_top();
            if (top == target) {
                distance = distFromSource[target];
                break;
            }

            G->forNeighborsOf(top, [&](node u, edgeweight w) {
                const double newDist = distFromSource[top] + w;
                if (newDist < distFromSource[u]) {
                    distFromSource[u] = newDist;
                    prio[u] = newDist + heu(u);
                    heap.update(u);
                    if (storePred) {
                        pred[u] = top;
                    }
                }
            });
        } while (!heap.empty());

        if (top != target) {
            WARN("Source cannot reach target!");
        } else if (storePred) {
            buildPath();
        }

        hasRun = true;
    }

private:
    // Priorities used for the heap
    std::vector<double> distFromSource, prio;
    // Lower bound of the distance to the target node
    Heuristic heu;

    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<double>> heap;
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_A_STAR_GENERAL_HPP_
