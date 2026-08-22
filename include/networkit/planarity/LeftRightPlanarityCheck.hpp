/*  LeftRightPlanarityCheck.hpp
 *
 *  Created on: 03.01.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_
#define NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_

#include <limits>
#include <stack>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

template <typename GraphT>
class GenericLeftRightPlanarityCheck final : public Algorithm {

public:
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

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
    GenericLeftRightPlanarityCheck(const GraphT &G);

    void run() override;

    bool isPlanar() const {
        assureFinished();
        return isGraphPlanar;
    }

private:
    count numberOfEdges;
    static constexpr count noneHeight{std::numeric_limits<count>::max()};
    static constexpr edgeid noneEdgeId{std::numeric_limits<edgeid>::max()};

    struct Interval {
        edgeid low{noneEdgeId};
        edgeid high{noneEdgeId};

        Interval() = default;
        Interval(edgeid lowId, edgeid highId) : low(lowId), high(highId) {}

        bool isEmpty() const { return low == noneEdgeId && high == noneEdgeId; }

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

    const ConflictPair NoneConflictPair{Interval{}, Interval{}};

    const GraphT *graph;
    bool isGraphPlanar{false};

    void dfsOrientation(NodeT startNode);
    bool dfsTesting(NodeT startNode);

    bool applyConstraints(edgeid edgeId, edgeid parentEdgeId);
    void removeBackEdges(edgeid edgeId, NodeT parentNode);
    void sortAdjacencyListByNestingDepth();
    bool conflicting(const Interval &interval, edgeid edgeId);
    count getLowestLowPoint(const ConflictPair &conflictPair);

    // DFS / lowpoint state
    std::vector<count> heights; // per node
    std::vector<NodeT> roots;   // roots of DFS forest

    std::vector<count> lowestPoint;
    std::vector<count> secondLowestPoint;
    std::vector<edgeid> ref;
    std::vector<edgeid> lowestPointEdge;
    std::vector<count> nestingDepth;
    std::vector<ConflictPair> stackBottom;

    // Per-node parent edge + parent node in DFS tree
    std::vector<edgeid> parentEdgeIds; // parent edge for each node (noneEdgeId if root)
    std::vector<NodeT> parentNodes;    // parent node for each node (none if root)

    // For each *underlying* edge ID, remember the DFS orientation's head (v)
    std::vector<NodeT> edgeEndpoints; // edgeEndpoints[eid] = v in oriented DFS edge (u -> v)

    // Conflict stack
    std::stack<ConflictPair> stack;
    // DFS graph with tree/back orientation
    GraphT dfsGraph;
};

using LeftRightPlanarityCheck = GenericLeftRightPlanarityCheck<Graph>;

} // namespace NetworKit

#include <networkit/planarity/LeftRightPlanarityCheckImpl.hpp>

#endif // NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_
