/*
 * DynDijkstra.hpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#ifndef NETWORKIT_DISTANCE_DYN_DIJKSTRA_HPP_
#define NETWORKIT_DISTANCE_DYN_DIJKSTRA_HPP_

#include <span>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/distance/DynSSSP.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Dynamic Dijkstra.
 */
class DynDijkstra final : public DynSSSP {

public:
    /**
     * Creates the object for @a G and source @a s.
     *
     * @param G The graph.
     * @param s The source node.
     * @param storePredecessors keep track of the lists of predecessors?
     */
    DynDijkstra(const Graph &G, node s, bool storePredecessors = true);

    // TODO the run method could take a vector of distances as an input and in that case just use
    // those distances instead of computing dijkstra from scratch
    void run() override;

    /** Updates the distances after an edge insertion.*/
    void update(GraphEvent e) override;

    /** Updates the distances after a batch of edge insertions.*/
    void updateBatch(std::span<const GraphEvent> batch) override;

private:
    enum Color { WHITE, GRAY, BLACK };
    std::vector<Color> color;
    static constexpr edgeweight infDist = std::numeric_limits<edgeweight>::max();
    static constexpr edgeweight distEpsilon = 0.000001;

    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<edgeweight>> heap;
    std::vector<edgeweight> updateDistances;
    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<edgeweight>> updateHeap;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_DYN_DIJKSTRA_HPP_
