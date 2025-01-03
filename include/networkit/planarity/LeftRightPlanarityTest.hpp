//
// Created by andreas on 03.01.25.
//

#ifndef LEFT_RIGHT_PLANARITY_TEST_HPP
#define LEFT_RIGHT_PLANARITY_TEST_HPP
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class LeftRightPlanarityTest final : public Algorithm {

public:
    LeftRightPlanarityTest(const Graph &graph): graph_(&graph){}
    void run() override;
	bool is_planar() const
    {
        return is_planar_;
    }

private:
    const Graph *graph_;
	bool is_planar_{};
    void dfsOrientation();
    void dfsTesting();
    bool applyConstraints();
    void removeBackEdges();
    std::vector<index> height;
    std::vector<count> lowPt;
    std::vector<count> lowPt2;
	std::unordered_map<Edge, Edge> lowPtEdge;
    std::unordered_map<Edge, count> nestingDepth;
    std::unordered_map<index, Edge> parentEdges;
    std::unordered_set<Edge> visitedEdges;
    Graph dfsGraph;

};
}
#endif //LEFT_RIGHT_PLANARITY_TEST_HPP
