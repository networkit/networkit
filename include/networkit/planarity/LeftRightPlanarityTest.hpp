//
// Created by andreas on 03.01.25.
//

#ifndef LEFT_RIGHT_PLANARITY_TEST_HPP
#define LEFT_RIGHT_PLANARITY_TEST_HPP
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

const Edge LRNoneEdge{none, none};

struct Interval
{
    Edge low{LRNoneEdge}; // Represents the lower bound of the interval
    Edge high{LRNoneEdge}; // Represents the upper bound of the interval

    // Default constructor
    Interval() = default;

    // Constructor with specific low and high values
    Interval(const Edge& low, const Edge& high)
        : low(low), high(high)
    {
    }

    bool is_empty() const
    {
        return low == LRNoneEdge && high == LRNoneEdge;
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
    std::vector<count> lowPoint;
    std::vector<count> secondLowestPoint;
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
