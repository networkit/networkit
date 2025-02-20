/*  FloydWarshall.hpp
 *
 *	Created on: 15.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_
#define NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
namespace NetworKit {

class FloydWarshall : public Algorithm {
public:
    FloydWarshall(const Graph &G);

    void run() override;

    edgeweight getDistance(node source, node target) const;

    const std::vector<std::vector<edgeweight>> & getAllDistances() const &;

    bool isNodeInNegativeCycle(node u) const;

    std::vector<node> getNodesOnShortestPath(node source, node target) const;

private:
    const Graph *graph;
    std::vector<std::vector<edgeweight>> distances;
    std::vector<uint8_t> nodesInNegativeCycle;
    std::vector<std::vector<node>> pathMatrix;
    void tagNegativeCycles();
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_
