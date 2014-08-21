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
	 * Get the index of node @a v in adjacency list of @a u in @a graph.
	 * @param graph The graph.
	 * @param u The start node of the edge.
	 * @param v The end node of the edge.
	 * @return The index of the edge.
	 */
	index getEdgeIdx(const Graph &graph, node u, node v) const;

	/**
	 * Performs a breadth-first search on @a graph from node @a source to find an augmenting path to @a sink respecting the @a flow.
	 * @param graph The graph.
	 * @param flow The current flow in the network.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param pred Used to store the path from @a source to @a sink.
	 * @return The gain in terms of flow.
	 */
	edgeweight BFS(const Graph &graph, std::vector<std::vector<edgeweight>> &flow, const node source, const node sink, std::vector<node> &pred) const;

	/**
	 * Performs the Edmonds-Karp algorithm on @a graph with @a source and @a sink. The flow is stored in @a flow.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param flow Used to store the flow on each edge.
	 * @return The maximum flow.
	 */
	edgeweight solveMaxFlow(const Graph &graph, const node source, const node sink, std::vector<std::vector<edgeweight>> &flow) const;

	/**
	 * Computes the source set and stores it in @a sourceSet.
	 * @param graph The graph.
	 * @param source The source node.
	 * @param sink The sink node.
	 * @param flow Used to store the flow on each edge.
	 * @param sourceSet Used to store the nodes belonging to the source set.
	 */
	void computeSourceSet(const Graph &graph, const node source, const node sink, const std::vector<std::vector<edgeweight>> &flow, std::vector<node> &sourceSet) const;

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
	edgeweight run(Graph &graph, const node source, const node sink, int &attribute_id) const;

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
	edgeweight run(Graph &graph, const node source, const node sink, std::vector<node> &sourceSet, int &attribute_id) const;


};

} /* namespace NetworKit */

#endif /* EDMONDSKARP_H_ */
