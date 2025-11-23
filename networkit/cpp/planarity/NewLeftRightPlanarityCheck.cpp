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

void NewLeftRightPlanarityCheck::run() {
    // Euler-criterion: non-planar if m > 3n - 6 for n > 2
    if (graph->numberOfNodes() > 2
        && graph->numberOfEdges() > 3 * graph->numberOfNodes() - 6) {
        hasRun = true;
        isGraphPlanar = false;
        return;
    }

    // Prepare per-node / per-edge state
    heights.assign(graph->upperNodeIdBound(), noneHeight);
    roots.clear();

    parentEdgeId.assign(graph->upperNodeIdBound(), noneEdgeId);
    parentNodes.assign(graph->upperNodeIdBound(), noneNode);

    edgeHead.assign(graph->upperEdgeIdBound(), noneNode);

    lowestPoint.clear();
    secondLowestPoint.clear();
    ref.clear();
    lowestPointEdge.clear();
    nestingDepth.clear();
    stackBottom.clear();
    edgeToNodesDEBUG = std::vector<std::pair<node, node>>(numberOfEdges);
    graph->forEdges([&](node u, node v, edgeweight w, edgeid id) {
        edgeToNodesDEBUG[id] = {u, v};
    });
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
        dfsGraph.sortNeighbors(currentNode,
                               [&](node neighbor1, node neighbor2) {
                                   const edgeid e1 = graph->edgeId(currentNode, neighbor1);
                                   const edgeid e2 = graph->edgeId(currentNode, neighbor2);

                                   auto it1 = nestingDepth.find(e1);
                                   auto it2 = nestingDepth.find(e2);

                                   if (it1 != nestingDepth.end() && it2 != nestingDepth.end()) {
                                       return it1->second < it2->second;
                                   }
                                   return false;
                               });
    });
}

bool NewLeftRightPlanarityCheck::conflicting(const Interval &interval, edgeid edgeId) {
    auto iteratorHigh = lowestPoint.find(interval.high);
    auto iteratorEdge = lowestPoint.find(edgeId);

    return !interval.isEmpty()
           && iteratorHigh != lowestPoint.end()
           && iteratorEdge != lowestPoint.end()
           && iteratorHigh->second > iteratorEdge->second;
}

bool NewLeftRightPlanarityCheck::applyConstraints(edgeid edgeId, edgeid parentEdgeId_) {
    ConflictPair tmpConflictPair{};

    // First phase: pop until stackBottom[edgeId], merging intervals on the right side.
    do {
        ConflictPair currentConflictPair = stack.top();
        stack.pop();

        if (!currentConflictPair.left.isEmpty()) {
            currentConflictPair.swap();
        }
        if (!currentConflictPair.left.isEmpty()) {
            // both sides non-empty -> non-planar
            return false;
        }

        auto rightLowIterator   = lowestPoint.find(currentConflictPair.right.low);
        auto parentEdgeIterator = lowestPoint.find(parentEdgeId_);

        if (rightLowIterator != lowestPoint.end()
            && parentEdgeIterator != lowestPoint.end()
            && rightLowIterator->second > parentEdgeIterator->second) {

            if (tmpConflictPair.right.isEmpty()) {
                tmpConflictPair.right = currentConflictPair.right;
            } else {
                ref[tmpConflictPair.right.low] = currentConflictPair.right.high;
            }

            tmpConflictPair.right.low = currentConflictPair.right.low;
        } else {
            // Attach to lowest point edge of the parent
            auto itLP = lowestPointEdge.find(parentEdgeId_);
            if (itLP != lowestPointEdge.end()) {
                ref[currentConflictPair.right.low] = itLP->second;
            }
        }
    } while (!stack.empty() && (stack.top() != stackBottom[edgeId]));

    // Second phase: pop while there are conflicts with edgeId, merging left/right intervals.
    while (!stack.empty()
           && (conflicting(stack.top().left, edgeId)
               || conflicting(stack.top().right, edgeId))) {

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
    return std::min(lowestPoint[conflictPair.right.low],
                    lowestPoint[conflictPair.left.low]);
}

void NewLeftRightPlanarityCheck::removeBackEdges(const edgeid edgeId, const node parentNode) {
    while (!stack.empty()
           && getLowestLowPoint(stack.top()) == heights[parentNode]) {
        stack.pop();
    }

    if (!stack.empty()) {
        auto conflictPair = stack.top();
        stack.pop();

        // Reduce left interval
        while (conflictPair.left.high != noneEdgeId
               && edgeHead[conflictPair.left.high] == parentNode) {
            auto it = ref.find(conflictPair.left.high);
            conflictPair.left.high = (it != ref.end()) ? it->second : noneEdgeId;
        }
        if (conflictPair.left.high == noneEdgeId
            && conflictPair.left.low != noneEdgeId) {
            ref[conflictPair.left.low] = conflictPair.right.low;
            conflictPair.left.low = noneEdgeId;
        }

        // Reduce right interval
        while (conflictPair.right.high != noneEdgeId
               && edgeHead[conflictPair.right.high] == parentNode) {
            auto it = ref.find(conflictPair.right.high);
            conflictPair.right.high = (it != ref.end()) ? it->second : noneEdgeId;
        }
        if (conflictPair.right.high == noneEdgeId
            && conflictPair.right.low != noneEdgeId) {
            ref[conflictPair.right.low] = conflictPair.left.low;
            conflictPair.right.low = noneEdgeId;
        }

        stack.push(conflictPair);
    }

    auto lpIt = lowestPoint.find(edgeId);
    if (!stack.empty()
        && lpIt != lowestPoint.end()
        && lpIt->second < heights[parentNode]) {

        const edgeid highestReturnEdgeLeft  = stack.top().left.high;
        const edgeid highestReturnEdgeRight = stack.top().right.high;

        edgeid chosen = noneEdgeId;

        if (highestReturnEdgeLeft != noneEdgeId) {
            if (highestReturnEdgeRight == noneEdgeId
                || lowestPoint[highestReturnEdgeLeft]
                       > lowestPoint[highestReturnEdgeRight]) {
                chosen = highestReturnEdgeLeft;
            } else {
                chosen = highestReturnEdgeRight;
            }
        } else {
            chosen = highestReturnEdgeRight;
        }

        if (chosen != noneEdgeId) {
            ref[edgeId] = chosen;
        }
    }
}

bool NewLeftRightPlanarityCheck::dfsTesting(node startNode) {
    std::stack<node> dfsStack;
    dfsStack.push(startNode);

    // Per-node neighbor iterators in DFS graph
    std::unordered_map<node, decltype(dfsGraph.neighborRange(startNode).begin())>
        neighborIterators;

    // We now deduplicate edges by edge IDs instead of Edge
    std::unordered_set<edgeid> preprocessedEdges;
    preprocessedEdges.reserve(numberOfEdges);

    auto processNeighborEdges =
        [&](node currentNode, bool &callRemoveBackEdges) -> bool {
            auto &neighborIterator = neighborIterators[currentNode];
            auto range             = dfsGraph.neighborRange(currentNode);

            while (neighborIterator != range.end()) {
                const node neighbor = *neighborIterator;
                const edgeid eid    = graph->edgeId(currentNode, neighbor);

                if (!preprocessedEdges.contains(eid)) {
                    stackBottom[eid] = stack.empty() ? NoneConflictPair : stack.top();

                    if (eid == parentEdgeId[neighbor]) {
                        // Tree edge: go deeper
                        dfsStack.push(currentNode);
                        dfsStack.push(neighbor);
                        preprocessedEdges.insert(eid);
                        callRemoveBackEdges = false;
                        ++neighborIterator;
                        return true; // more work later
                    }

                    lowestPointEdge[eid] = eid;
                    stack.emplace(Interval{}, Interval(eid, eid));
                    preprocessedEdges.insert(eid);
                }

                auto currentEdgeIterator = lowestPoint.find(eid);
                if (currentEdgeIterator != lowestPoint.end()
                    && currentEdgeIterator->second < heights[currentNode]) {

                    if (neighbor == *range.begin()) {
                        // First neighbor: propagate lowestPointEdge to parent edge
                        lowestPointEdge[parentEdgeId[currentNode]] =
                            lowestPointEdge[eid];
                    } else if (!applyConstraints(eid, parentEdgeId[currentNode])) {
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

        const edgeid parentEid = parentEdgeId[currentNode];
        bool callRemoveBackEdges{true};

        if (!neighborIterators.contains(currentNode)) {
            neighborIterators[currentNode] =
                dfsGraph.neighborRange(currentNode).begin();
        }

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

        const edgeid parentEid = parentEdgeId[currentNode];

        for (node neighbor : graph->neighborRange(currentNode)) {
            const edgeid eid = graph->edgeId(currentNode, neighbor);

            if (!preprocessedEdges.contains(eid)) {
                if (dfsGraph.hasEdge(currentNode, neighbor)
                    || dfsGraph.hasEdge(neighbor, currentNode)) {
                    continue;
                }

                dfsGraph.addEdge(currentNode, neighbor);
                edgeHead[eid] = neighbor;

                lowestPoint[eid]       = heights[currentNode];
                secondLowestPoint[eid] = heights[currentNode];

                if (heights[neighbor] == noneHeight) {
                    // Tree edge
                    parentEdgeId[neighbor] = eid;
                    parentNodes[neighbor]   = currentNode;

                    heights[neighbor] = heights[currentNode] + 1;

                    dfsStack.push(currentNode);
                    dfsStack.push(neighbor);

                    preprocessedEdges.insert(eid);
                    break;
                }

                // Back edge: lowpoint is ancestor's height
                lowestPoint[eid] = heights[neighbor];
            }

            nestingDepth[eid] = 2 * lowestPoint[eid];
            if (secondLowestPoint[eid] < heights[currentNode]) {
                nestingDepth[eid] += 1;
            }

            if (parentEid != noneEdgeId) {
                if (lowestPoint[eid] < lowestPoint[parentEid]) {
                    secondLowestPoint[parentEid] =
                        std::min(lowestPoint[parentEid],
                                 secondLowestPoint[eid]);
                    lowestPoint[parentEid] = lowestPoint[eid];
                } else if (lowestPoint[eid] > lowestPoint[parentEid]) {
                    secondLowestPoint[parentEid] =
                        std::min(secondLowestPoint[parentEid],
                                 lowestPoint[eid]);
                } else {
                    secondLowestPoint[parentEid] =
                        std::min(secondLowestPoint[parentEid],
                                 secondLowestPoint[eid]);
                }
            }
        }
    } while (!dfsStack.empty());
}

} // namespace NetworKit
//
// Created by andreas on 22.11.25.
//
