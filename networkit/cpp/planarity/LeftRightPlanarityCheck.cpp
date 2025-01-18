/*  LeftRightPlanarityCheck.cpp
 *
 *	Created on: 03.01.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/planarity/LeftRightPlanarityCheck.hpp>

namespace NetworKit {

void LeftRightPlanarityCheck::run() {
    // Euler-criterion
    if (graph_->numberOfNodes() > 2 && graph_->numberOfEdges() > 3 * graph_->numberOfNodes() - 6) {
        hasRun = true;
        isPlanar_ = false;
        return;
    }

    heights.assign(graph_->upperNodeIdBound(), noneHeight);
    graph_->forNodes([&](const node currentNode) {
        if (heights[currentNode] == noneHeight) {
            heights[currentNode] = 0;
            roots.push_back(currentNode);
            this->dfsOrientation(currentNode);
        }
    });

    sortAdjacencyListByNestingDepth();
    isPlanar_ =
        std::ranges::all_of(roots, [this](const auto rootNode) { return dfsTesting(rootNode); });
    hasRun = true;
}

void LeftRightPlanarityCheck::sortAdjacencyListByNestingDepth() {

    dfsGraph.forNodes([&](const node currentNode) {
        dfsGraph.sortNeighbors(currentNode, [&](const node neighbor1, const node neighbor2) {
            if (auto it1 = nestingDepth.find(Edge(currentNode, neighbor1)),
                it2 = nestingDepth.find(Edge(currentNode, neighbor2));
                it1 != nestingDepth.end() && it2 != nestingDepth.end()) {
                return it1->second < it2->second;
            }
            return false;
        });
    });
}

bool LeftRightPlanarityCheck::conflicting(const Interval &interval, const Edge &edge) {
    return !interval.isEmpty() && lowestPoint.contains(interval.high) && lowestPoint.contains(edge)
           && lowestPoint[interval.high] > lowestPoint[edge];
}

bool LeftRightPlanarityCheck::applyConstraints(const Edge edge, const Edge parentEdge) {
    auto tmpConflictPair = ConflictPair{};
    do {
        auto currentConflictPair = stack.top();
        stack.pop();
        if (!currentConflictPair.left.isEmpty()) {
            currentConflictPair.swap();
        }
        if (!currentConflictPair.left.isEmpty()) {
            return false;
        }
        if (lowestPoint.contains(currentConflictPair.right.low) && lowestPoint.contains(parentEdge)
            && lowestPoint[currentConflictPair.right.low] > lowestPoint[parentEdge]) {
            if (tmpConflictPair.right.isEmpty()) {
                tmpConflictPair.right = currentConflictPair.right;
            } else {
                ref[tmpConflictPair.right.low] = currentConflictPair.right.high;
            }

            tmpConflictPair.right.low = currentConflictPair.right.low;
        } else {
            ref[currentConflictPair.right.low] = lowestPointEdge[parentEdge];
        }
    } while (!stack.empty() && (stack.top() != stackBottom[edge]));

    while (!stack.empty()
           && (conflicting(stack.top().left, edge) || conflicting(stack.top().right, edge))) {
        auto currentConflictPair = stack.top();
        stack.pop();
        if (conflicting(currentConflictPair.right, edge)) {
            currentConflictPair.swap();
        }
        if (conflicting(currentConflictPair.right, edge)) {
            return false;
        }
        ref[tmpConflictPair.right.low] = currentConflictPair.right.high;
        if (currentConflictPair.right.low != noneEdge) {
            tmpConflictPair.right = currentConflictPair.right;
        }
        if (tmpConflictPair.left.isEmpty()) {
            tmpConflictPair.left = currentConflictPair.left;
        } else {
            ref[tmpConflictPair.left.low] = currentConflictPair.left.high;
        }
        tmpConflictPair.left.low = currentConflictPair.left.low;
    }

    if (!tmpConflictPair.left.isEmpty() || !tmpConflictPair.right.isEmpty())
        stack.push(tmpConflictPair);
    return true;
}

count LeftRightPlanarityCheck::getLowestLowPoint(const ConflictPair &conflictPair) {
    if (conflictPair.left.isEmpty())
        return lowestPoint[conflictPair.right.low];
    if (conflictPair.right.isEmpty())
        return lowestPoint[conflictPair.left.low];
    return std::min(lowestPoint[conflictPair.right.low], lowestPoint[conflictPair.left.low]);
}

void LeftRightPlanarityCheck::removeBackEdges(Edge edge) {
    auto parentNode = edge.u;
    while (!stack.empty() && getLowestLowPoint(stack.top()) == heights[parentNode]) {
        stack.pop();
    }

    if (!stack.empty()) {
        auto conflictPair = stack.top();
        stack.pop();
        while (conflictPair.left.high != noneEdge && conflictPair.left.high.v == parentNode) {
            auto it = ref.find(conflictPair.left.high);
            conflictPair.left.high = (it != ref.end()) ? it->second : noneEdge;
        }
        if (conflictPair.left.high == noneEdge && conflictPair.left.low != noneEdge) {
            ref[conflictPair.left.low] = conflictPair.right.low;
            conflictPair.left.low = noneEdge;
        }
        while (conflictPair.right.high != noneEdge && conflictPair.right.high.v == parentNode) {
            auto it = ref.find(conflictPair.right.high);
            conflictPair.right.high = (it != ref.end()) ? it->second : noneEdge;
        }

        if (conflictPair.right.high == noneEdge && conflictPair.right.low != noneEdge) {
            ref[conflictPair.right.low] = conflictPair.left.low;
            conflictPair.right.low = noneEdge;
        }
        stack.push(conflictPair);
    }

    if (!stack.empty() && lowestPoint[edge] < heights[parentNode]) {
        auto highestReturnEdgeLeft = stack.top().left.high;
        auto highestReturnEdgeRight = stack.top().right.high;
        if (highestReturnEdgeLeft != noneEdge
            && (highestReturnEdgeRight != noneEdge
                || lowestPoint[highestReturnEdgeLeft] > lowestPoint[highestReturnEdgeRight])) {
            ref[edge] = highestReturnEdgeLeft;
        } else {
            ref[edge] = highestReturnEdgeRight;
        }
    }
}

bool LeftRightPlanarityCheck::dfsTesting(node startNode) {
    std::stack<node> dfsStack;
    dfsStack.emplace(startNode);
    auto neighborIterators =
        std::unordered_map<node, decltype(dfsGraph.neighborRange(startNode).begin())>{};

    auto preprocessedEdges = std::unordered_set<Edge>{};
    while (!dfsStack.empty()) {
        const auto currentNode = dfsStack.top();
        dfsStack.pop();
        const auto parentEdge = parentEdges[currentNode];
        bool callRemoveBackEdges{true};
        if (auto it = neighborIterators.find(currentNode); it == neighborIterators.end()) {
            neighborIterators[currentNode] = dfsGraph.neighborRange(currentNode).begin();
        }
        auto &neighborIterator = neighborIterators[currentNode];
        while (neighborIterator != dfsGraph.neighborRange(currentNode).end()) {
            const auto neighbor = *neighborIterator;
            auto currentEdge = Edge(currentNode, neighbor);
            if (!preprocessedEdges.contains(currentEdge)) {
                stackBottom[currentEdge] = stack.empty() ? NoneConflictPair : stack.top();
                if (currentEdge == parentEdges[neighbor]) {
                    dfsStack.emplace(currentNode);
                    dfsStack.emplace(neighbor);
                    preprocessedEdges.insert(currentEdge);
                    callRemoveBackEdges = false;
                    break;
                }
                lowestPointEdge[currentEdge] = currentEdge;
                stack.emplace(Interval{}, Interval(currentEdge, currentEdge));
            }

            if (lowestPoint.contains(currentEdge)
                && lowestPoint[currentEdge] < heights[currentNode]) {
                if (neighbor == *dfsGraph.neighborRange(currentNode).begin()) {
                    lowestPointEdge[parentEdge] = lowestPointEdge[currentEdge];
                } else {
                    if (!applyConstraints(currentEdge, parentEdge))
                        return false;
                }
            }
            ++neighborIterator;
        }

        if (callRemoveBackEdges) {
            if (parentEdge != noneEdge)
                removeBackEdges(parentEdge);
        }
    }
    return true;
}

void LeftRightPlanarityCheck::dfsOrientation(const node startNode) {

    std::stack<node> dfsStack;
    dfsStack.emplace(startNode);
    auto preprocessedEdges = std::unordered_set<Edge>{};
    while (!dfsStack.empty()) {
        const auto currentNode = dfsStack.top();
        dfsStack.pop();
        const auto parentEdge = parentEdges[currentNode];
        for (const auto neighbor : graph_->neighborRange(currentNode)) {

            const auto currentEdge = Edge(currentNode, neighbor);
            if (!preprocessedEdges.contains(currentEdge)) {
                if (dfsGraph.hasEdge(currentNode, neighbor)
                    || dfsGraph.hasEdge(neighbor, currentNode))
                    continue;
                dfsGraph.addEdge(currentNode, neighbor);
                lowestPoint[currentEdge] = heights[currentNode];
                secondLowestPoint[currentEdge] = heights[currentNode];
                if (heights[neighbor] == noneHeight) // Tree edge
                {
                    parentEdges[neighbor] = currentEdge;
                    heights[neighbor] = heights[currentNode] + 1;
                    dfsStack.emplace(currentNode);
                    dfsStack.emplace(neighbor);
                    preprocessedEdges.insert(currentEdge);
                    break;
                }
                // back edge
                lowestPoint[currentEdge] = heights[neighbor];
            }

            nestingDepth[currentEdge] = 2 * lowestPoint[currentEdge];
            if (secondLowestPoint[currentEdge] < heights[currentNode]) {
                nestingDepth[currentEdge] += 1;
            }

            if (parentEdge != noneEdge) {
                if (lowestPoint[currentEdge] < lowestPoint[parentEdge]) {
                    secondLowestPoint[parentEdge] =
                        std::min(lowestPoint[parentEdge], secondLowestPoint[currentEdge]);
                    lowestPoint[parentEdge] = lowestPoint[currentEdge];
                } else if (lowestPoint[currentEdge] > lowestPoint[parentEdge]) {
                    secondLowestPoint[parentEdge] =
                        std::min(secondLowestPoint[parentEdge], lowestPoint[currentEdge]);
                } else {
                    secondLowestPoint[parentEdge] =
                        std::min(secondLowestPoint[parentEdge], secondLowestPoint[currentEdge]);
                }
            }
        }
    }
}

} // namespace NetworKit
