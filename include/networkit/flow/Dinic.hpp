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
        Dinic(const Graph &G, node s, node t): graph(&G), source(s), target(t) {
            residual = Graph(graph->upperNodeIdBound(), true, true);
        }
        void run() override;
    private:
        void buildResidual();
        bool bfs();
        edgeweight dfs(node u, edgeweight flow);
        node source;
        node target;
        const Graph *graph;
        edgeweight maxFlow{};
        Graph residual;
        std::vector<int> level;
    std::vector<int> ptr;

};
}
#endif //NETWORKIT_FLOW_DINIC_HPP_
