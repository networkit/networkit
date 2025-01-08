//
// Created by andreas on 03.01.25.
//
#include <networkit/planarity/LeftRightPlanarityTest.hpp>

namespace NetworKit {

void LeftRightPlanarityTest::run() {
    // Euler-criterion
    if (graph_->numberOfNodes() > 2 && graph_->numberOfEdges() > 3 * graph_->numberOfNodes() - 6) {
        is_planar_ = false;
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
    is_planar_ = true;
    for (const auto rootNode : roots) {
        if (!dfsTesting(rootNode)) {
            is_planar_ = false;
            break;
        }
    }
}

void LeftRightPlanarityTest::sortAdjacencyListByNestingDepth() {

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

bool LeftRightPlanarityTest::conflicting(const Interval &interval, const Edge &edge) {
    return !interval.is_empty() && lowestPoint.contains(interval.high) && lowestPoint.contains(edge) && lowestPoint[interval.high] > lowestPoint[edge];
}

bool LeftRightPlanarityTest::applyConstraints(const Edge edge, const Edge parentEdge) {
    auto help_conflict_pair = ConflictPair{};
    do {
        auto current_conflict_pair = stack.top();
        stack.pop();
        if (!current_conflict_pair.left.is_empty()) {
            current_conflict_pair.swap();
        }
        if (!current_conflict_pair.left.is_empty()) {
            return false;
        }
        if (lowestPoint.contains(current_conflict_pair.right.low)
            && lowestPoint.contains(parentEdge)
            && lowestPoint[current_conflict_pair.right.low] > lowestPoint[parentEdge]) {
            if (help_conflict_pair.right.is_empty()) {
                help_conflict_pair.right = current_conflict_pair.right;
            } else {
                ref[help_conflict_pair.right.low] = current_conflict_pair.right.high;
            }

            help_conflict_pair.right.low = current_conflict_pair.right.low;
        } else {
            ref[current_conflict_pair.right.low] = lowestPointEdge[parentEdge];
        }
    } while (!stack.empty() && (stack.top() != stackBottom[edge]));

    while (!stack.empty()
           && (conflicting(stack.top().left, edge) || conflicting(stack.top().right, edge))) {
        auto current_conflict_pair = stack.top();
        stack.pop();
        if (conflicting(current_conflict_pair.right, edge)) {
            current_conflict_pair.swap();
        }
        if (conflicting(current_conflict_pair.right, edge)) {
            return false;
        }
        ref[help_conflict_pair.right.low] = current_conflict_pair.right.high;
        if (current_conflict_pair.right.low != noneEdge) {
            help_conflict_pair.right = current_conflict_pair.right;
        }
        if (help_conflict_pair.left.is_empty()) {
            help_conflict_pair.left = current_conflict_pair.left;
        } else {
            ref[help_conflict_pair.left.low] = current_conflict_pair.left.high;
        }
        help_conflict_pair.left.low = current_conflict_pair.left.low;
    }

    if (!help_conflict_pair.left.is_empty() || !help_conflict_pair.right.is_empty())
        stack.push(help_conflict_pair);
    return true;
}

count LeftRightPlanarityTest::getLowestLowPoint(const ConflictPair &conflict_pair) {
    if (conflict_pair.left.is_empty())
        return lowestPoint[conflict_pair.right.low];
    if (conflict_pair.right.is_empty())
        return lowestPoint[conflict_pair.left.low];
    return std::min(lowestPoint[conflict_pair.right.low], lowestPoint[conflict_pair.left.low]);
}

void LeftRightPlanarityTest::removeBackEdges(Edge edge) {
    auto parent_node = edge.u;
    while (!stack.empty() && getLowestLowPoint(stack.top()) == heights[parent_node]) {
        stack.pop();
    }

    if (!stack.empty()) {
        auto conflict_pair = stack.top();
        stack.pop();
        while (conflict_pair.left.high != noneEdge && conflict_pair.left.high.v == parent_node) {
            auto it = ref.find(conflict_pair.left.high);
            conflict_pair.left.high = (it != ref.end()) ? it->second : noneEdge;
        }
        if (conflict_pair.left.high == noneEdge && conflict_pair.left.low != noneEdge) {
            ref[conflict_pair.left.low] = conflict_pair.right.low;
            conflict_pair.left.low = noneEdge;
        }
        while (conflict_pair.right.high != noneEdge && conflict_pair.right.high.v == parent_node) {
            auto it = ref.find(conflict_pair.right.high);
            conflict_pair.right.high = (it != ref.end()) ? it->second : noneEdge;
        }

        if (conflict_pair.right.high == noneEdge && conflict_pair.right.low != noneEdge) {
            ref[conflict_pair.right.low] = conflict_pair.left.low;
            conflict_pair.right.low = noneEdge;
        }
        stack.push(conflict_pair);
    }

    if (!stack.empty() && lowestPoint[edge] < heights[parent_node]) {
        auto highest_return_edge_left = stack.top().left.high;
        auto highest_return_edge_right = stack.top().right.high;
        if (highest_return_edge_left != noneEdge
            && (highest_return_edge_right != noneEdge
                || lowestPoint[highest_return_edge_left]
                       > lowestPoint[highest_return_edge_right])) {
            ref[edge] = highest_return_edge_left;
        } else {
            ref[edge] = highest_return_edge_right;
        }
    }
}

bool LeftRightPlanarityTest::dfsTesting(node startNode) {
    std::stack<node> dfs_stack;
    dfs_stack.emplace(startNode);
    auto neighborIterators = std::unordered_map<node, decltype(dfsGraph.neighborRange(startNode).begin())>{};

    auto preprocessed_edges = std::unordered_set<Edge>{};
    while (!dfs_stack.empty()) {
        const auto currentNode = dfs_stack.top();
        dfs_stack.pop();
        const auto parentEdge = parentEdges[currentNode];
        bool callRemoveBackEdges{true};
        if (auto it = neighborIterators.find(currentNode); it == neighborIterators.end()) {
            neighborIterators[currentNode] = dfsGraph.neighborRange(currentNode).begin();
        }
        auto& neighborIterator = neighborIterators[currentNode];
        while (neighborIterator != dfsGraph.neighborRange(currentNode).end()) {
            const auto neighbor = *neighborIterator;
            auto currentEdge = Edge(currentNode, neighbor);
            if (!preprocessed_edges.contains(currentEdge)) {
                stackBottom[currentEdge] = stack.empty() ? NoneConflictPair : stack.top();
                if (currentEdge == parentEdges[neighbor]) {
                    dfs_stack.emplace(currentNode);
                    dfs_stack.emplace(neighbor);
                    preprocessed_edges.insert(currentEdge);
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

void LeftRightPlanarityTest::dfsOrientation(const node startNode) {

    std::stack<node> dfs_stack;
    dfs_stack.emplace(startNode);
    auto preprocessed_edges = std::unordered_set<Edge>{};
    while (!dfs_stack.empty()) {
        const auto currentNode = dfs_stack.top();
        dfs_stack.pop();
        const auto parentEdge = parentEdges[currentNode];
        for (const auto neighbor : graph_->neighborRange(currentNode)) {

            const auto currentEdge = Edge(currentNode, neighbor);
            if (!preprocessed_edges.contains(currentEdge)) {
                if (dfsGraph.hasEdge(currentNode, neighbor) || dfsGraph.hasEdge(neighbor, currentNode))
                    continue;
                dfsGraph.addEdge(currentNode, neighbor);
                lowestPoint[currentEdge] = heights[currentNode];
                secondLowestPoint[currentEdge] = heights[currentNode];
                if (heights[neighbor] == noneHeight) // Tree edge
                {
                    parentEdges[neighbor] = currentEdge;
                    heights[neighbor] = heights[currentNode] + 1;
                    dfs_stack.emplace(currentNode);
                    dfs_stack.emplace(neighbor);
                    preprocessed_edges.insert(currentEdge);
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
