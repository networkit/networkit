//
// Created by andreas on 22.11.25.
//

#ifndef NETWORKIT_NEWLEFTRIGHTPLANARITYCHECK_HPP
#define NETWORKIT_NEWLEFTRIGHTPLANARITYCHECK_HPP
/*  NewLeftRightPlanarityCheck.hpp
 *
 *  Created on: 03.01.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_
#define NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_

#include <limits>
#include <stack>
#include <unordered_map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class NewLeftRightPlanarityCheck final : public Algorithm {

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
    NewLeftRightPlanarityCheck(const Graph &G) : graph(&G) {
        if (G.isDirected()) {
            throw std::runtime_error("The graph is not an undirected graph.");
        }
        if (!G.hasEdgeIds()) {
            throw std::runtime_error("The graph has no edge ids..");
        }


        numberOfEdges = graph->numberOfEdges();

        lowestPoint.reserve(numberOfEdges);
        secondLowestPoint.reserve(numberOfEdges);
        ref.reserve(numberOfEdges);
        lowestPointEdge.reserve(numberOfEdges);
        nestingDepth.reserve(numberOfEdges);
        stackBottom.reserve(numberOfEdges);

        // dfsGraph: directed view of DFS tree + back edges
        dfsGraph = Graph(graph->numberOfNodes(), /*weighted=*/false,
                         /*directed=*/true, /*edgesIndexed=*/false);
    }

    void run() override;

    bool isPlanar() const {
        assureFinished();
        return isGraphPlanar;
    }

private:
    // We still keep a none-Edge as a *value* sentinel (but never as a key).
    static const Edge noneEdge;

    count numberOfEdges;
    static constexpr count noneHeight{std::numeric_limits<count>::max()};
    static constexpr edgeid noneEdgeId{std::numeric_limits<edgeid>::max()};
    static constexpr node noneNode{std::numeric_limits<node>::max()};

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

    const Graph *graph;
    bool isGraphPlanar{false};

    // Algorithm phases
    void dfsOrientation(node startNode);
    bool dfsTesting(node startNode);

    // All these now work with edge IDs
    bool applyConstraints(const edgeid edgeId, const edgeid parentEdgeId);
    void removeBackEdges(const edgeid edgeId, const node parentNode);
    void sortAdjacencyListByNestingDepth();
    bool conflicting(const Interval &interval, edgeid edgeId);
    count getLowestLowPoint(const ConflictPair &conflictPair);

    // DFS / lowpoint state
    std::vector<count> heights;        // per node
    std::vector<node> roots;           // roots of DFS forest

    // All maps keyed by edgeid instead of Edge now
    std::unordered_map<edgeid, count> lowestPoint;
    std::unordered_map<edgeid, count> secondLowestPoint;
    std::unordered_map<edgeid, edgeid> ref;
    std::unordered_map<edgeid, edgeid> lowestPointEdge;
    std::vector<count> nestingDepth;
    std::unordered_map<edgeid, ConflictPair> stackBottom;

    // Per-node parent edge + parent node in DFS tree (by ID, not Edge)
    std::vector<edgeid> parentEdgeId;  // parent edge for each node (noneEdgeId if root)
    std::vector<node> parentNodes;      // parent node for each node (noneNode if root)

    // For each *underlying* edge ID, remember the DFS orientation's head (v)
    std::vector<node> edgeHead;        // edgeHead[eid] = v in oriented DFS edge (u -> v)

    // Conflict stack
    std::stack<ConflictPair> stack;
    std::vector<std::pair<node, node>> edgeToNodesDEBUG;
    // DFS graph with tree/back orientation
    Graph dfsGraph;
};

} // namespace NetworKit

#endif // NETWORKIT_PLANARITY_LEFT_RIGHT_PLANARITY_CHECK_HPP_

#endif // NETWORKIT_NEWLEFTRIGHTPLANARITYCHECK_HPP
