/*  LeftRightPlanarityCheck.hpp
 *
 *	Created on: 03.01.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_
#define NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class LeftRightPlanarityCheck final : public Algorithm {

public:
    /**
     * Implements the left-right planarity test as described in
     * [citation](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7963e9feffe1c9362eb1a69010a5139d1da3661e).
     * This algorithm determines whether a graph is planar, i.e., whether it can be drawn on a plane
     * without any edges crossing. For an overview of planar graphs, refer to:
     * https://en.wikipedia.org/wiki/Planar_graph
     *
     * The algorithm achieves (almost) linear runtime complexity. The only non-linear component
     * arises from sorting the nodes of the depth-first search tree.
     *
     * @param G The input graph to test for planarity. The graph should be undirected.
     * @throws std::runtime_error if graph is not an undirected graph
     */
    LeftRightPlanarityCheck(const Graph &G) : graph(&G) {
        if (G.isDirected()) {
            throw std::runtime_error("The graph is not an undirected graph.");
        }
        dfsGraph = Graph(graph->numberOfNodes(), false, true, false);
    }

    /**
     * Executes the left-right planarity test on the input graph.
     * This method performs all necessary computations to determine
     * whether the graph is planar and prepares the result for retrieval
     * via the `isPlanar()` method.
     */
    void run() override;

    /**
     * Returns whether the input graph is planar.
     * The result is only valid after the `run()` method has been called.
     *
     * @return True if the graph is planar, false otherwise.
     * @throws std::runtime_error if called before `run()` has been executed.
     */
    bool isPlanar() const {
        assureFinished();
        return isGraphPlanar;
    }

private:
    static const Edge<node> noneEdge;
    static constexpr count noneHeight{std::numeric_limits<count>::max()};

    struct Interval {
        Edge<node> low{noneEdge};
        Edge<node> high{noneEdge};

        Interval() : low{noneEdge}, high{noneEdge} {};
        Interval(const Edge<node> &low, const Edge<node> &high) : low(low), high(high) {}
        bool isEmpty() const { return low == noneEdge && high == noneEdge; }

        friend bool operator==(const Interval &lhs, const Interval &rhs) {
            return lhs.low == rhs.low && lhs.high == rhs.high;
        }
    };

    struct ConflictPair {
        Interval left{};
        Interval right{};

        ConflictPair() = default;
        ConflictPair(const Interval &left, const Interval &right) : left(left), right(right) {}

        void swap() { std::swap(left, right); }

        friend bool operator==(const ConflictPair &lhs, const ConflictPair &rhs) {
            return lhs.left == rhs.left && lhs.right == rhs.right;
        }
    };
    const ConflictPair NoneConflictPair{Interval(), Interval()};

    const Graph *graph;
    bool isGraphPlanar{};
    void dfsOrientation(node startNode);
    bool dfsTesting(node startNode);
    bool applyConstraints(const Edge<node> &edge, const Edge<node> &parentEdge);
    void removeBackEdges(const Edge<node> &edge);
    void sortAdjacencyListByNestingDepth();
    bool conflicting(const Interval &interval, const Edge<node> &edge);
    count getLowestLowPoint(const ConflictPair &conflictPair);
    std::vector<count> heights;
    std::unordered_map<Edge<node>, count> lowestPoint;
    std::unordered_map<Edge<node>, count> secondLowestPoint;
    std::unordered_map<Edge<node>, Edge<node>> ref;
    std::vector<node> roots;
    std::unordered_map<Edge<node>, Edge<node>> lowestPointEdge;
    std::unordered_map<Edge<node>, count> nestingDepth;
    std::unordered_map<index, Edge<node>> parentEdges;
    std::stack<ConflictPair> stack;
    std::unordered_map<Edge<node>, ConflictPair> stackBottom;
    Graph dfsGraph;
};
} // namespace NetworKit
#endif // NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_
