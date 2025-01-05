//
// Created by andreas on 03.01.25.
//

#ifndef LEFT_RIGHT_PLANARITY_TEST_HPP
#define LEFT_RIGHT_PLANARITY_TEST_HPP
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

const Edge noneEdge{none, none};

struct Interval
{
    Edge low{noneEdge}; // Represents the lower bound of the interval
    Edge high{noneEdge}; // Represents the upper bound of the interval

    // Default constructor
    Interval() = default;

    // Constructor with specific low and high values
    Interval(const Edge& low, const Edge& high)
        : low(low), high(high)
    {
    }

    bool is_empty() const
    {
        return low == noneEdge && high == noneEdge;
    }
};

inline bool operator==(const Interval& lhs, const Interval& rhs)
{
    return lhs.low == rhs.low && lhs.high == rhs.high;
}


struct ConflictPair
{
    Interval left{}; // Left interval of edges
    Interval right{}; // Right interval of edges

    ConflictPair() = default;

    // Constructor with initial intervals
    ConflictPair(const Interval& left, const Interval& right)
        : left(left), right(right)
    {
    }

    void swap()
    {
        std::swap(left, right);
    }
};

const ConflictPair NoneConflictPair{Interval{}, Interval{}};

inline  bool operator==(const ConflictPair& lhs, const ConflictPair& rhs)
{
    return lhs.left == rhs.left && lhs.right == rhs.right;
}


class LeftRightPlanarityTest final : public Algorithm {

public:
    LeftRightPlanarityTest(const Graph &graph): graph_(&graph){}

    void run() override;

    void initialization();

    bool is_planar() const
    {
        return is_planar_;
    }



private:
    constexpr count noneHeight{std::numeric_limits<count>::max()};
    const Graph *graph_;
	bool is_planar_{};
    void dfsOrientation(node startNode);
    bool dfsTesting(node startNode);
    bool applyConstraints(Edge edge, Edge parentEdge);
    void removeBackEdges(Edge edge);
    std::vector<count> heights;
    std::unordered_map<Edge, count> lowestPoint;
    std::unordered_map<Edge, count> secondLowestPoint;
    std::vector<node> roots;
	std::unordered_map<Edge, Edge> lowPointEdge;
    std::unordered_map<Edge, count> nestingDepth;
    std::unordered_map<index, Edge> parentEdges;
    std::unordered_set<Edge> visitedEdges;
    std::stack<ConflictPair> stack;
    std::unordered_map<Edge, ConflictPair> stackBottom;
    Graph dfsGraph;

};
}
#endif //LEFT_RIGHT_PLANARITY_TEST_HPP
