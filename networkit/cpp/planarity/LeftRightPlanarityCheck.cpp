/*  LeftRightPlanarityCheck.cpp
 *
 *	Created on: 03.01.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/planarity/LeftRightPlanarityCheck.hpp>

namespace NetworKit {

const Edge LeftRightPlanarityCheck::noneEdge{};

void LeftRightPlanarityCheck::run() {
    // Euler-criterion
    if (graph->numberOfNodes() > 2 && graph->numberOfEdges() > 3 * graph->numberOfNodes() - 6) {
        hasRun = true;
        isGraphPlanar = false;
        return;
    }

    heights.assign(graph->upperNodeIdBound(), noneHeight);
    graph->forNodes([&](node currentNode) {
        if (heights[currentNode] == noneHeight) {
            heights[currentNode] = 0;
            roots.push_back(currentNode);
            this->dfsOrientation(currentNode);
        }
    });

    sortAdjacencyListByNestingDepth();
    isGraphPlanar =
        std::ranges::all_of(roots, [this](node rootNode) { return dfsTesting(rootNode); });
    hasRun = true;
}

void LeftRightPlanarityCheck::sortAdjacencyListByNestingDepth() {
    dfsGraph.forNodes([&](node currentNode) {
        dfsGraph.sortNeighbors(currentNode, [&](node neighbor1, node neighbor2) {
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
    auto iteratorHigh = lowestPoint.find(interval.high);
    auto iteratorEdge = lowestPoint.find(edge);
    return !interval.isEmpty() && iteratorHigh != lowestPoint.end()
           && iteratorEdge != lowestPoint.end() && iteratorHigh->second > iteratorEdge->second;
}

bool LeftRightPlanarityCheck::applyConstraints(const Edge &edge, const Edge &parentEdge) {
    ConflictPair tmpConflictPair{};
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
        auto parentEdgeIterator = lowestPoint.find(parentEdge);
        if (rightLowIterator != lowestPoint.end() && parentEdgeIterator != lowestPoint.end()
            && rightLowIterator->second > parentEdgeIterator->second) {
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

void LeftRightPlanarityCheck::removeBackEdges(const Edge &edge) {
    const node parentNode = edge.u;
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
        const Edge highestReturnEdgeLeft = stack.top().left.high;
        const Edge highestReturnEdgeRight = stack.top().right.high;
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
    dfsStack.push(startNode);
    std::unordered_map<node, decltype(dfsGraph.neighborRange(startNode).begin())> neighborIterators;
    std::unordered_set<Edge> preprocessedEdges;

    auto processNeighborEdges = [&](node currentNode, bool &callRemoveBackEdges) -> bool {
        auto &neighborIterator = neighborIterators[currentNode];
        while (neighborIterator != dfsGraph.neighborRange(currentNode).end()) {
            const node neighbor = *neighborIterator;
            const Edge currentEdge(currentNode, neighbor);

            if (!preprocessedEdges.contains(currentEdge)) {
                stackBottom[currentEdge] = stack.empty() ? NoneConflictPair : stack.top();
                if (currentEdge == parentEdges[neighbor]) {
                    dfsStack.push(currentNode);
                    dfsStack.push(neighbor);
                    preprocessedEdges.insert(currentEdge);
                    callRemoveBackEdges = false;
                    return true; // Indicate further processing needed
                }

                lowestPointEdge[currentEdge] = currentEdge;
                stack.emplace(Interval{}, Interval(currentEdge, currentEdge));
            }

            if (auto currentEdgeIterator = lowestPoint.find(currentEdge);
                currentEdgeIterator != lowestPoint.end()
                && currentEdgeIterator->second < heights[currentNode]) {

                if (neighbor == *dfsGraph.neighborRange(currentNode).begin()) {
                    lowestPointEdge[parentEdges[currentNode]] = lowestPointEdge[currentEdge];
                } else if (!applyConstraints(currentEdge, parentEdges[currentNode])) {
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
        const Edge parentEdge = parentEdges[currentNode];
        bool callRemoveBackEdges{true};

        if (auto it = neighborIterators.find(currentNode); it == neighborIterators.end()) {
            neighborIterators[currentNode] = dfsGraph.neighborRange(currentNode).begin();
        }

        if (!processNeighborEdges(currentNode, callRemoveBackEdges)) {
            return false;
        }

        if (callRemoveBackEdges && parentEdge != noneEdge) {
            removeBackEdges(parentEdge);
        }
    } while (!dfsStack.empty());

    return true;
}

void LeftRightPlanarityCheck::dfsOrientation(node startNode) {
    std::stack<node> dfsStack;
    dfsStack.push(startNode);
    std::unordered_set<Edge> preprocessedEdges;
    do {
        const node currentNode = dfsStack.top();
        dfsStack.pop();
        const Edge parentEdge = parentEdges[currentNode];
        for (node neighbor : graph->neighborRange(currentNode)) {

            const Edge currentEdge = Edge(currentNode, neighbor);
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
                    dfsStack.push(currentNode);
                    dfsStack.push(neighbor);
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
    } while (!dfsStack.empty());
}

} // namespace NetworKit
