#ifndef NETWORKIT_CPP_DISTANCE_VERTEX_DIAMETER_BOUND_HPP_
#define NETWORKIT_CPP_DISTANCE_VERTEX_DIAMETER_BOUND_HPP_

#include <limits>
#include <span>

#include <networkit/Globals.hpp>
#include <networkit/graph/GraphConcepts.hpp>

namespace NetworKit::DiameterDetails {

/**
 * Upper bound on the vertex diameter of the connected component of a node in an undirected graph,
 * given the distances from that node: a shortest path between two other nodes is at most as long
 * as the path through the node, so the two largest finite distances add up to at least the
 * diameter.
 */
template <GraphEdgeWeight EdgeWeightT>
count vertexDiameterBound(std::span<const EdgeWeightT> distances) {
    EdgeWeightT largest = 0;
    EdgeWeightT secondLargest = 0;
    for (EdgeWeightT dist : distances) {
        if (dist == std::numeric_limits<EdgeWeightT>::max()) {
            continue;
        }
        if (dist >= largest) {
            secondLargest = largest;
            largest = dist;
        } else if (dist > secondLargest) {
            secondLargest = dist;
        }
    }
    // Count the nodes, not the edges.
    return static_cast<count>(largest + secondLargest) + 1;
}

} // namespace NetworKit::DiameterDetails

#endif // NETWORKIT_CPP_DISTANCE_VERTEX_DIAMETER_BOUND_HPP_
