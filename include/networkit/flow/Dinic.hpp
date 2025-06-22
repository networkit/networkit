/*  Dinic.hpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_FLOW_DINIC_HPP_
#define NETWORKIT_FLOW_DINIC_HPP_
#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class Dinic final : public Algorithm {

public:
    Dinic(const Graph &G, node s, node t);
    void run() override;

private:
    void buildResidual();
    bool get_parents_bfs();
    edgeweight blocking_path();
    node source;
    node target;
    const Graph *graph;
    edgeweight maxFlow{};
    Graph residualGraph;
    std::vector<std::deque<node>> parents;
    std::vector<int> levels;
};
} // namespace NetworKit
#endif // NETWORKIT_FLOW_DINIC_HPP_
