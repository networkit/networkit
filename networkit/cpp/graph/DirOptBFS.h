/*
 * DirOptBFS.h
 *
 *  Created on: Feb 21, 2015
 *      Author: Maximilian Vogel
 */

#ifndef DIROPTBFS_H_
#define DIROPTBFS_H_

#include "Graph.h"
#include "SSSP.h"

namespace NetworKit {

/**
 * @ingroup graph
 * The BFS class is used to do a breadth-first search on a Graph from a given source node.
 * Idea of the direction optimizing BFS came from: http://www.cs.berkeley.edu/~sbeamer/beamer-sc2012.pdf
 */
class DirOptBFS : public SSSP {

friend class DynBFS;

private:
	std::vector<bool> frontier; //TODO: doc
	std::vector<bool> next; // TODO: doc

	/**
	 * Indicates if the BFS is finished yet.
	 */
	bool hasQueuedNodes;

	/**
	 * Indicates the direction of the next step.
	 */
	bool topdown;

	/**
	 * Tuning parameter needed when determining whether it's necessary to switch from top-down to bottom-up.
	 */
	count alpha;

	/**
	 * Tuning parameter needed when determining whether it's necessary to switch from bottom-up to top-down.
	 */
	count beta;

	/**
	 * Number of edges to be visited from the frontier (to be computed on the fly).
	 */
	count m_f;

	/**
	 * Number of edges to be visited from unvisited nodes (should decrease monotonously?).
	 */
	count m_u;

	/**
	 * Track number of nodes in the frontier.
	 */
	count n_f;

	/**
	 * The other part of the equation with the nodes.
	 */
	count rhs_C_BT;

public:
	/**
	 * Constructs the BFS class for @a G and source node @a source.
	 *
	 * @param G The graph.
	 * @param source The source node of the breadth-first search.
	 * @param storePaths	store paths and number of paths?
	 * @param storeStack	maintain a stack of nodes in decreasing order of distance
	 */
	DirOptBFS(const Graph& G, node source, bool storePaths=true, bool storeStack=false, count alphe = 12, count beta = 24, node target = none);

	/**
	 * Breadth-first search from @a source.
	 * @return Vector of unweighted distances from node @a source, i.e. the
	 * length (number of edges) of the shortest path from @a source to any other node.
	 */
	virtual void run() override;

};

} /* namespace NetworKit */
#endif /* BFS_H_ */
