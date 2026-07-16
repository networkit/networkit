/*
 * SSSPImpl.hpp
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#ifndef NETWORKIT_DISTANCE_SSSP_IMPL_HPP_
#define NETWORKIT_DISTANCE_SSSP_IMPL_HPP_

#include <stack>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

template <class GraphT>
SingleSourceShortestPaths<GraphT>::SingleSourceShortestPaths(const GraphT &G, NodeT source,
                                                             bool storePaths,
                                                             bool storeNodesSortedByDistance,
                                                             NodeT target)
    : Algorithm(), G(&G), source(source), target(target), storePaths(storePaths),
      storeNodesSortedByDistance(storeNodesSortedByDistance) {}

template <class GraphT>
auto SingleSourceShortestPaths<GraphT>::getPath(NodeT t, bool forward) const -> std::vector<NodeT> {
    if (!storePaths) {
        throw std::runtime_error("paths have not been stored");
    }
    std::vector<NodeT> path;
    if (previous[t].empty()) { // t is not reachable from source
        WARN("there is no path from ", source, " to ", t);
        return path;
    }
    NodeT v = t;
    while (v != source) {
        path.push_back(v);
        v = previous[v].front();
    }
    path.push_back(source);

    if (forward) {
        std::reverse(path.begin(), path.end());
    }
    return path;
}

template <class GraphT>
auto SingleSourceShortestPaths<GraphT>::getPaths(NodeT t, bool forward) const
    -> std::set<std::vector<NodeT>> {

    std::set<std::vector<NodeT>> paths;
    if (previous[t].empty()) { // t is not reachable from source
        WARN("there is no path from ", source, " to ", t);
        return paths;
    }

    std::vector<NodeT> targetPredecessors = previous[t];

#pragma omp parallel for schedule(dynamic)
    for (omp_index i = 0; i < static_cast<omp_index>(targetPredecessors.size()); ++i) {

        std::stack<std::vector<NodeT>> stack;
        std::vector<std::vector<NodeT>> currPaths;

        NodeT pred = targetPredecessors[i];
        if (pred == source) {
            currPaths.push_back({t, pred});
        } else {
            stack.push({t, pred});
        }

        while (!stack.empty()) {

            NodeT topPathLastNode = stack.top().back();

            if (topPathLastNode == source) {
                currPaths.push_back(stack.top());
                stack.pop();
                continue;
            }

            std::vector<NodeT> topPath = stack.top();
            stack.pop();

            std::vector<NodeT> currPredecessors = previous[topPath.back()];

            for (NodeT currPredecessor : currPredecessors) {
                std::vector<NodeT> suffix(topPath);
                suffix.push_back(currPredecessor);
                stack.push(suffix);
            }
        }

#pragma omp critical
        paths.insert(currPaths.begin(), currPaths.end());
    }

    if (forward) {
        std::set<std::vector<NodeT>> reversedPaths;
        for (std::vector<NodeT> path : paths) {
            std::reverse(std::begin(path), std::end(path));
            reversedPaths.insert(path);
        }
        paths = reversedPaths;
    }

    return paths;
}

template <class GraphT>
auto SingleSourceShortestPaths<GraphT>::getNodesSortedByDistance() const
    -> const std::vector<NodeT> & {
    if (!storeNodesSortedByDistance) {
        throw std::runtime_error("Nodes sorted by distance have not been stored. Set "
                                 "storeNodesSortedByDistance in the constructor to true to enable "
                                 "this behaviour.");
    }

    assert(!nodesSortedByDistance.empty());
    return nodesSortedByDistance;
}

#endif // NETWORKIT_DISTANCE_SSSP_IMPL_HPP_

} /* namespace NetworKit */
