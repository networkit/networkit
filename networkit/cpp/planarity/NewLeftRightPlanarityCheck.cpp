/*  NewLeftRightPlanarityCheck.cpp
 *
 *  Created on: 03.01.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <algorithm>
#include <unordered_set>

#include <networkit/planarity/NewLeftRightPlanarityCheck.hpp>

namespace NetworKit {

const Edge NewLeftRightPlanarityCheck::noneEdge{};

NewLeftRightPlanarityCheck::NewLeftRightPlanarityCheck(const Graph &G) : graph(&G) {
    if (G.isDirected()) {
        throw std::runtime_error("The graph is not an undirected graph.");
    }

    numberOfEdges = graph->numberOfEdges();
    nestingDepth.resize(numberOfEdges, -1);
    lowestPointEdge.resize(numberOfEdges, noneEdgeId);

    lowestPoint.reserve(numberOfEdges);
    secondLowestPoint.resize(numberOfEdges, none);
    ref.reserve(numberOfEdges);

    stackBottom.reserve(numberOfEdges);

    // dfsGraph: directed view of DFS tree + back edges
    dfsGraph = Graph(graph->numberOfNodes(), /*weighted=*/false,
                     /*directed=*/true, /*edgesIndexed=*/false);
}

void NewLeftRightPlanarityCheck::run() {
    // Euler-criterion: non-planar if m > 3n - 6 for n > 2
    if (graph->numberOfNodes() > 2 && graph->numberOfEdges() > 3 * graph->numberOfNodes() - 6) {
        hasRun = true;
        isGraphPlanar = false;
        return;
    }
    if (!graph->hasEdgeIds()) {
        // Logical constness: indexing does not change topology, only adds IDs.
        const_cast<Graph *>(graph)->indexEdges();
    }
    // Prepare per-node / per-edge state
    heights.assign(graph->upperNodeIdBound(), noneHeight);
    roots.clear();

    parentEdgeIds.assign(graph->upperNodeIdBound(), noneEdgeId);
    parentNodes.assign(graph->upperNodeIdBound(), noneNode);

    edgeEndpoints.assign(graph->upperEdgeIdBound(), noneNode);

    lowestPoint.clear();
    ref.clear();
    stackBottom.clear();
    edgeToNodesDEBUG = std::vector<std::pair<node, node>>(numberOfEdges);
    graph->forEdges(
        [&](node u, node v, edgeweight w, edgeid id) { edgeToNodesDEBUG[id] = {u, v}; });
    // DFS orientation
    graph->forNodes([&](node currentNode) {
        if (heights[currentNode] == noneHeight) {
            heights[currentNode] = 0;
            roots.push_back(currentNode);
            this->dfsOrientation(currentNode);
        }
    });

    // Sort DFS adjacency by nesting depth
    sortAdjacencyListByNestingDepth();

    // Planarity testing DFS
    isGraphPlanar =
        std::ranges::all_of(roots, [this](node rootNode) { return dfsTesting(rootNode); });

    hasRun = true;
}

void NewLeftRightPlanarityCheck::sortAdjacencyListByNestingDepth() {
    dfsGraph.forNodes([&](node currentNode) {
        dfsGraph.sortNeighbors(currentNode, [&](node neighbor1, node neighbor2) {
            const edgeid e1 = graph->edgeId(currentNode, neighbor1);
            const edgeid e2 = graph->edgeId(currentNode, neighbor2);
            return nestingDepth[e1] < nestingDepth[e2];
        });
    });
}

bool NewLeftRightPlanarityCheck::conflicting(const Interval &interval, edgeid edgeId) {
    auto iteratorHigh = lowestPoint.find(interval.high);
    auto iteratorEdge = lowestPoint.find(edgeId);

    return !interval.isEmpty() && iteratorHigh != lowestPoint.end()
           && iteratorEdge != lowestPoint.end() && iteratorHigh->second > iteratorEdge->second;
}

bool NewLeftRightPlanarityCheck::applyConstraints(edgeid edgeId, edgeid parentEdgeId) {
    ConflictPair tmpConflictPair{};

    // First phase: pop until stackBottom[edgeId], merging intervals on the right side.
    do {
        ConflictPair currentConflictPair = stack.top();
        stack.pop();

        if (!currentConflictPair.left.isEmpty()) {
            currentConflictPair.swap();
        }
        if (!currentConflictPair.left.isEmpty()) {
            return false;
        }

        auto rightLowIterator = lowestPoint.find(currentConflictPair.right.low);
        auto parentEdgeIterator = lowestPoint.find(parentEdgeId);

        if (rightLowIterator != lowestPoint.end() && parentEdgeIterator != lowestPoint.end()
            && rightLowIterator->second > parentEdgeIterator->second) {

            if (tmpConflictPair.right.isEmpty()) {
                tmpConflictPair.right = currentConflictPair.right;
            } else {
                ref[tmpConflictPair.right.low] = currentConflictPair.right.high;
            }

            tmpConflictPair.right.low = currentConflictPair.right.low;
        } else {
            ref[currentConflictPair.right.low] = lowestPointEdge[parentEdgeId];
        }
    } while (!stack.empty() && (stack.top() != stackBottom[edgeId]));

    // Second phase: pop while there are conflicts with edgeId, merging left/right intervals.
    while (!stack.empty()
           && (conflicting(stack.top().left, edgeId) || conflicting(stack.top().right, edgeId))) {

        auto currentConflictPair = stack.top();
        stack.pop();

        if (conflicting(currentConflictPair.right, edgeId)) {
            currentConflictPair.swap();
        }
        if (conflicting(currentConflictPair.right, edgeId)) {
            // Still conflicting after swap -> non-planar
            return false;
        }

        // NOTE: this is intentionally exactly like the original code:
        // it writes to ref[tmpConflictPair.right.low] even if tmpConflictPair.right was empty.
        ref[tmpConflictPair.right.low] = currentConflictPair.right.high;

        if (currentConflictPair.right.low != noneEdgeId) {
            tmpConflictPair.right = currentConflictPair.right;
        }

        if (tmpConflictPair.left.isEmpty()) {
            tmpConflictPair.left = currentConflictPair.left;
        } else {
            ref[tmpConflictPair.left.low] = currentConflictPair.left.high;
        }
        tmpConflictPair.left.low = currentConflictPair.left.low;
    }

    if (!tmpConflictPair.left.isEmpty() || !tmpConflictPair.right.isEmpty()) {
        stack.push(tmpConflictPair);
    }

    return true;
}

count NewLeftRightPlanarityCheck::getLowestLowPoint(const ConflictPair &conflictPair) {
    if (conflictPair.left.isEmpty()) {
        return lowestPoint[conflictPair.right.low];
    }
    if (conflictPair.right.isEmpty()) {
        return lowestPoint[conflictPair.left.low];
    }
    return std::min(lowestPoint[conflictPair.right.low], lowestPoint[conflictPair.left.low]);
}

void NewLeftRightPlanarityCheck::removeBackEdges(const edgeid edgeId, const node parentNode) {
    while (!stack.empty() && getLowestLowPoint(stack.top()) == heights[parentNode]) {
        stack.pop();
    }

    if (!stack.empty()) {
        auto conflictPair = stack.top();
        stack.pop();

        // Reduce left interval
        while (conflictPair.left.high != noneEdgeId
               && edgeEndpoints[conflictPair.left.high] == parentNode) {
            auto it = ref.find(conflictPair.left.high);
            conflictPair.left.high = (it != ref.end()) ? it->second : noneEdgeId;
        }
        if (conflictPair.left.high == noneEdgeId && conflictPair.left.low != noneEdgeId) {
            ref[conflictPair.left.low] = conflictPair.right.low;
            conflictPair.left.low = noneEdgeId;
        }

        // Reduce right interval
        while (conflictPair.right.high != noneEdgeId
               && edgeEndpoints[conflictPair.right.high] == parentNode) {
            auto it = ref.find(conflictPair.right.high);
            conflictPair.right.high = (it != ref.end()) ? it->second : noneEdgeId;
        }
        if (conflictPair.right.high == noneEdgeId && conflictPair.right.low != noneEdgeId) {
            ref[conflictPair.right.low] = conflictPair.left.low;
            conflictPair.right.low = noneEdgeId;
        }

        stack.push(conflictPair);
    }

    if (!stack.empty() && lowestPoint[edgeId] < heights[parentNode]) {
        const edgeid highestReturnEdgeLeft = stack.top().left.high;
        const edgeid highestReturnEdgeRight = stack.top().right.high;

        if (highestReturnEdgeLeft != noneEdgeId
            && (highestReturnEdgeRight != noneEdgeId
                || lowestPoint[highestReturnEdgeLeft] > lowestPoint[highestReturnEdgeRight])) {
            ref[edgeId] = highestReturnEdgeLeft;
        } else {
            ref[edgeId] = highestReturnEdgeRight;
        }
    }
}

bool NewLeftRightPlanarityCheck::dfsTesting(node startNode) {
    std::stack<node> dfsStack;
    dfsStack.push(startNode);

    // Per-node neighbor iterators in DFS graph
    std::unordered_map<node, decltype(dfsGraph.neighborRange(startNode).begin())> neighborIterators;

    // We now deduplicate edges by edge IDs instead of Edge
    std::unordered_set<edgeid> preprocessedEdges;
    preprocessedEdges.reserve(numberOfEdges);

    auto processNeighborEdges = [&](node currentNode, bool &callRemoveBackEdges) -> bool {
        auto &neighborIterator = neighborIterators[currentNode];
        while (neighborIterator != dfsGraph.neighborRange(currentNode).end()) {
            const node neighbor = *neighborIterator;
            const edgeid currentEdgeId = graph->edgeId(currentNode, neighbor);

            if (!preprocessedEdges.contains(currentEdgeId)) {
                stackBottom[currentEdgeId] = stack.empty() ? NoneConflictPair : stack.top();

                if (currentEdgeId == parentEdgeIds[neighbor]) {
                    // Tree edge: go deeper
                    dfsStack.push(currentNode);
                    dfsStack.push(neighbor);
                    preprocessedEdges.insert(currentEdgeId);
                    callRemoveBackEdges = false;
                    return true; // more work later
                }

                lowestPointEdge[currentEdgeId] = currentEdgeId;
                stack.emplace(Interval{}, Interval(currentEdgeId, currentEdgeId));
            }

            auto currentEdgeIterator = lowestPoint.find(currentEdgeId);
            if (currentEdgeIterator != lowestPoint.end()
                && currentEdgeIterator->second < heights[currentNode]) {

                if (neighbor == *dfsGraph.neighborRange(currentNode).begin()) {
                    // First neighbor: propagate lowestPointEdge to parent edge
                    lowestPointEdge[parentEdgeIds[currentNode]] = lowestPointEdge[currentEdgeId];
                } else if (!applyConstraints(currentEdgeId, parentEdgeIds[currentNode])) {
                    return false;
                }
            }

            ++neighborIterator;
        }

        return true;
    };

    // Main DFS loop
    do {
        const node currentNode = dfsStack.top();
        dfsStack.pop();

        const edgeid parentEid = parentEdgeIds[currentNode];
        bool callRemoveBackEdges{true};

        if (auto it = neighborIterators.find(currentNode); it == neighborIterators.end()) {
            neighborIterators[currentNode] = dfsGraph.neighborRange(currentNode).begin();
        }
        // analysis when currentNode is 4 the 2nd time in 'testNonPlanarCompleteBipartiteGraphK3_3', 9 stackBottom are correct, 3 stack entries are correct
        // 6 parentEdgeIds are correct: lowestPointEdge has different length
        if (!processNeighborEdges(currentNode, callRemoveBackEdges)) {
            return false;
        }

        if (callRemoveBackEdges && parentEid != noneEdgeId) {
            removeBackEdges(parentEid, parentNodes[currentNode]);
        }

    } while (!dfsStack.empty());

    return true;
}

void NewLeftRightPlanarityCheck::dfsOrientation(node startNode) {
    std::stack<node> dfsStack;
    dfsStack.push(startNode);

    // Deduplicate per *underlying* edge, not per oriented Edge struct
    std::unordered_set<edgeid> preprocessedEdges;
    preprocessedEdges.reserve(numberOfEdges);

    do {
        const node currentNode = dfsStack.top();
        dfsStack.pop();

        const edgeid parentEdgeId = parentEdgeIds[currentNode];

        for (node neighbor : graph->neighborRange(currentNode)) {
            const edgeid edgeId = graph->edgeId(currentNode, neighbor);

            if (!preprocessedEdges.contains(edgeId)) {
                if (dfsGraph.hasEdge(currentNode, neighbor)
                    || dfsGraph.hasEdge(neighbor, currentNode)) {
                    continue;
                }

                dfsGraph.addEdge(currentNode, neighbor);
                edgeEndpoints[edgeId] = neighbor;

                lowestPoint[edgeId] = heights[currentNode];
                secondLowestPoint[edgeId] = heights[currentNode];

                if (heights[neighbor] == noneHeight) {
                    // Tree edge
                    parentEdgeIds[neighbor] = edgeId;
                    parentNodes[neighbor] = currentNode;

                    heights[neighbor] = heights[currentNode] + 1;

                    dfsStack.push(currentNode);
                    dfsStack.push(neighbor);

                    preprocessedEdges.insert(edgeId);
                    break;
                }

                // Back edge: lowpoint is ancestor's height
                lowestPoint[edgeId] = heights[neighbor];
            }

            nestingDepth[edgeId] = 2 * lowestPoint[edgeId];
            if (secondLowestPoint[edgeId] < heights[currentNode]) {
                nestingDepth[edgeId] += 1;
            }

            if (parentEdgeId != noneEdgeId) {
                if (lowestPoint[edgeId] < lowestPoint[parentEdgeId]) {
                    secondLowestPoint[parentEdgeId] =
                        std::min(lowestPoint[parentEdgeId], secondLowestPoint[edgeId]);
                    lowestPoint[parentEdgeId] = lowestPoint[edgeId];
                } else if (lowestPoint[edgeId] > lowestPoint[parentEdgeId]) {
                    secondLowestPoint[parentEdgeId] =
                        std::min(secondLowestPoint[parentEdgeId], lowestPoint[edgeId]);
                } else {
                    secondLowestPoint[parentEdgeId] =
                        std::min(secondLowestPoint[parentEdgeId], secondLowestPoint[edgeId]);
                }
            }
        }
    } while (!dfsStack.empty());
}

} // namespace NetworKit
//
// Created by andreas on 22.11.25.
//
