//
// Created by andreas on 03.01.25.
//
#include <networkit/planarity/LeftRightPlanarityTest.hpp>

namespace NetworKit {

void LeftRightPlanarityTest::initialization() {
    const auto n = graph_->upperNodeIdBound();
    heights.assign(n, noneHeight);
}

void LeftRightPlanarityTest::run() {
    // Euler-critirium

    initialization();
    if (graph_->numberOfNodes() > 2 && graph_->numberOfEdges() > 3 * graph_->numberOfNodes() - 6) {
        is_planar_ = false;
        return;
    }
    graph_->forNodes([&](const node currentNode) {
        if (heights[currentNode] == noneHeight) {
            heights[currentNode] = 0;
            roots.push_back(currentNode);
            this->dfsOrientation(currentNode);
        }
    });

    sortAdjacencyListByNestingDepth();
    is_planar_ = true;
    for(const auto rootNode : roots) {
        if(!dfsTesting(rootNode))
            is_planar_ = false;
    }
}

void LeftRightPlanarityTest::sortAdjacencyListByNestingDepth() {

    dfsGraph.forNodes([&](const node currentNode) {
        std::sort(dfsGraph.neighborRange(currentNode).begin(),
                  dfsGraph.neighborRange(currentNode).end(), [&](const node a, const node b) {
                      return nestingDepth.at(Edge(currentNode, a))
                             < nestingDepth.at(Edge(currentNode, b));
                  });
    });
}

bool LeftRightPlanarityTest::dfsTesting(node startNode) {
    std::stack<std::pair<node, Graph::NeighborIterator>> dfs_stack;
    dfs_stack.emplace(startNode, graph_->neighborRange(startNode).begin());

    auto preprocessed_edges = std::unordered_set<Edge>{};
    while (!dfs_stack.empty()) {
        const auto currentNode = dfs_stack.top().first;
        auto &neighborIterator = dfs_stack.top().second;
        dfs_stack.pop();
        const auto parentEdge = parentEdges[currentNode];
        bool callRemoveBackEdges{true};
        auto processed_neighbors = std::unordered_set<node>{};

        for (; neighborIterator != graph_->neighborRange(currentNode).end(); ++neighborIterator) {
            const auto neighbor = *neighborIterator;
            if (processed_neighbors.contains(neighbor))
                continue;
            processed_neighbors.insert(neighbor);

            auto currentEdge = Edge(currentNode, neighbor);
            if (!preprocessed_edges.contains(currentEdge)) {
                stackBottom[currentEdge] = stack.empty() ? NoneConflictPair : stack.top();
                if (currentEdge == parentEdges[neighbor]) {
                    dfs_stack.emplace(currentNode, graph_->neighborRange(currentNode).begin());
                    dfs_stack.emplace(neighbor, graph_->neighborRange(neighbor).begin());
                    preprocessed_edges.insert(currentEdge);
                    callRemoveBackEdges = false;
                    break;
                }
                lowPointEdge[currentEdge] = currentEdge;
                stack.emplace(Interval{}, Interval(currentEdge, currentEdge));
            }

            if (lowestPoint[currentEdge] < heights[currentNode]) {
                if (neighbor == *dfsGraph.neighborRange(currentNode).begin()) {
                    lowPointEdge[parentEdge] = lowPointEdge[currentEdge];
                } else {
                    if (!applyConstraints(currentEdge, parentEdge))
                        return false;
                }
            }
        }

        if (callRemoveBackEdges) {
            if (parentEdge != noneEdge)
                removeBackEdges(parentEdge);
        }
    }
    return true;
}

void LeftRightPlanarityTest::dfsOrientation(const node startNode) {

    std::stack<std::pair<node, Graph::NeighborIterator>> dfs_stack;
    dfs_stack.emplace(startNode, graph_->neighborRange(startNode).begin());
    auto preprocessed_edges = std::unordered_set<Edge>{};

    while (!dfs_stack.empty()) {
        const auto currentNode = dfs_stack.top().first;
        auto &neighborIterator = dfs_stack.top().second;
        dfs_stack.pop();
        const auto parentEdge = parentEdges[currentNode];
        auto processedNeighbors = std::unordered_set<node>{};

        for (; neighborIterator != graph_->neighborRange(currentNode).end(); ++neighborIterator) {
            const auto neighbor = *neighborIterator;
            if (processedNeighbors.contains(neighbor))
                continue;
            processedNeighbors.insert(neighbor);
            const auto currentEdge = Edge(currentNode, neighbor);
            if (!preprocessed_edges.contains(currentEdge)) {
                const auto currentReversedEdge = Edge(neighbor, currentNode);
                if (visitedEdges.contains(currentEdge)
                    || visitedEdges.contains(currentReversedEdge))
                    continue;
                visitedEdges.insert(currentEdge);
                dfsGraph.addEdge(currentNode, neighbor);
                lowestPoint[currentEdge] = heights[currentNode];
                secondLowestPoint[currentEdge] = heights[currentNode];
                if (heights[neighbor] == noneHeight) // Tree edge
                {
                    parentEdges[neighbor] = currentEdge;
                    heights[neighbor] = heights[currentNode] + 1;
                    dfs_stack.emplace(currentNode, graph_->neighborRange(currentNode).begin());
                    dfs_stack.emplace(neighbor, graph_->neighborRange(neighbor).begin());
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
