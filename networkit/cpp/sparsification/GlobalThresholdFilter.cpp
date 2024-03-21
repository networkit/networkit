/*
 * GlobalThresholdFilter.cpp
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/sparsification/GlobalThresholdFilter.hpp>

namespace NetworKit {

GlobalThresholdFilter::GlobalThresholdFilter(const Graph &graph,
                                             const std::vector<double> &attribute, double threshold,
                                             bool above)
    : graph(&graph), attribute(attribute), threshold(threshold), above(above) {}

Graph GlobalThresholdFilter::calculate() {
    if (!graph->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }
    if (!graph->isDirected())
        return calculateUndirected();

    // Create an edge-less graph.
    GraphBuilder builder(graph->upperNodeIdBound(), graph->isWeighted(), graph->isDirected());

    // Re-add the edges of the sparsified graph.
    graph->balancedParallelForNodes([&](const node u) {
        // add each edge in both directions
        graph->forEdgesOf(u, [&](const node u, const node v, const edgeweight ew,
                                 const edgeid eid) {
            if ((above && attribute[eid] >= threshold) || (!above && attribute[eid] <= threshold)) {
                builder.addHalfEdge(u, v, ew);
            }
        });
    });

    auto sGraph = builder.completeGraph();
    // WARNING: removeNode() must not be called in parallel (writes on vector<bool> and does
    // non-atomic decrement of number of nodes)!
    sGraph.forNodes([&](const node u) {
        if (!graph->hasNode(u)) {
            sGraph.removeNode(u);
        }
    });

    return sGraph;
}

Graph GlobalThresholdFilter::calculateUndirected() {
    // Create an edge-less graph.
    Graph sGraph(graph->upperNodeIdBound(), graph->isWeighted(), graph->isDirected());
    count edgeCount = 0;
    count selfLoops = 0;

    // Re-add the edges of the sparsified graph.
    graph->balancedParallelForNodes([&](const node u) {
        // add each edge in both directions
        graph->forEdgesOf(u, [&](const node u, const node v, const edgeweight ew,
                                 const edgeid eid) {
            if ((above && attribute[eid] >= threshold) || (!above && attribute[eid] <= threshold)) {
                sGraph.addPartialEdge(unsafe, u, v, ew);
#pragma omp atomic
                edgeCount++;
#pragma omp atomic
                selfLoops += (u == v);
            }
        });
    });
    edgeCount = (edgeCount + selfLoops) / 2;
    sGraph.setEdgeCount(unsafe, edgeCount);
    sGraph.setNumberOfSelfLoops(unsafe, selfLoops);
    // WARNING: removeNode() must not be called in parallel (writes on vector<bool> and does
    // non-atomic decrement of number of nodes)!
    sGraph.forNodes([&](const node u) {
        if (!graph->hasNode(u)) {
            sGraph.removeNode(u);
        }
    });

    return sGraph;
}

} /* namespace NetworKit */
