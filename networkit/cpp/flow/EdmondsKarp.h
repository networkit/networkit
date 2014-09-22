/*
 * EdmondsKarp.h
 *
 *  Created on: 11.06.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef EDMONDSKARP_H_
#define EDMONDSKARP_H_

#include "../graph/Graph.h"
#include "../algebraic/Matrix.h"
#include <vector>

namespace NetworKit {

/**
 * The EdmondsKarp class implements the maximum flow algorithm by Edmonds and Karp.
 */
class EdmondsKarp {
private:
	/**
	 * Performs a breadth-first search on @a graph from node @a source to find an augmenting path to @a sink respecting the @a flow.
	 * @param graph The graph.
	 * @param flow The current flow in the network.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param pred Used to store the path from @a source to @a sink.
	 * @return The gain in terms of flow.
	 */
	edgeweight BFS(const NetworKit::Graph &graph, std::vector< NetworKit::edgeweight > &flow, std::vector< NetworKit::edgeweight > &residFlow, node source, node sink, std::vector< NetworKit::node > &pred) const;

	/**
	 * Performs the Edmonds-Karp algorithm on @a graph with @a source and @a sink. The flow is stored in @a flow.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param flow Used to store the flow on each edge.
	 * @return The maximum flow.
	 */
	edgeweight solveMaxFlow(const NetworKit::Graph &graph, const node source, const node sink, std::vector< edgeweight > &flow) const;

	/**
	 * Computes the source set and stores it in @a sourceSet.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param flow Used to store the flow on each edge.
	 * @param sourceSet Used to store the nodes belonging to the source set.
	 */
	void computeSourceSet(const Graph &graph, const node source, const node sink, const std::vector<edgeweight> &flow, std::vector<node> &sourceSet) const;

public:
	/**
	 * Returns the maximum flow on @a graph with given @a source and @a sink.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @return The maximum flow.
	 */
	edgeweight run(const Graph &graph, const node source, const node sink) const;

	/**
	 * Returns the maximum flow on @a graph with given @a source and @a sink and stores the nodes belonging to the source set in @a sourceSet.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.#
	 * @param sourceSet Used to store the nodes belonging to the source set.
	 * @return The maximum flow.
	 */
	edgeweight run(const Graph &graph, const node source, const node sink, std::vector<node> &sourceSet) const;

	/**
	 * Returns the maximum flow on @a graph with given @a source and @a sink and stores the final flow as edge attributes with @a attribute_id.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param attribute_id The attribute id for getting the flow as edge attribute.
	 * @return The maximum flow.
	 */
	edgeweight run(const Graph &graph, const node source, const node sink, std::vector< edgeweight > &flow) const;

	/**
	 * Returns the maximum flow on @a graph with given @a source and @a sink. The nodes belonging to the source set are stored in @a sourceSet
	 * and the final flow is stored as edge attributes with @a attribute_id.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param sourceSet Used to store the nodes belonging to the source set.
	 * @param attribute_id The attribute id for getting the flow as edge attribute.
	 * @return The maximum flow.
	 */
	edgeweight run(const Graph &graph, node source, node sink, std::vector< node > &sourceSet, std::vector< edgeweight > &flow) const;


};

} /* namespace NetworKit */

#endif /* EDMONDSKARP_H_ */
