/*
 * ParallelConnectedComponents.h
 *
 *  Created on: Sep 17, 2013
 *      Author: Maximilian Vogel
 */

#ifndef PARALLELCONNECTEDCOMPONENTS_H_
#define PARALLELCONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../graph/BFS.h"

namespace NetworKit {

/**
 * Determines the connected components of an undirected graph.
 */
class ParallelConnectedComponents {
public:
	ParallelConnectedComponents();
	virtual ~ParallelConnectedComponents();
	/**
	 * This method determines the connected components for the graph g.
	 *
	 * @param[in]	G	graph 
	 */
	void run(const Graph& G);

	/**
	 * This method returns the number of connected components.
	 */
	count numberOfComponents();

	/**
	 * This method returns the the component in which node query is situated.
	 *
	 * @param[in]	query	the node whose component is asked for
	 */
	count componentOfNode(node u);


	/** 
	 * Return a vector which contains the
	 * component label indexed by each node
	 */
	std::vector<index> getComponentData();

	/**
	 * This method returns a component consisting of a vector of nodes
	 *
	 * @param[in]	component	the index of the component that is asked for
	 */
	std::vector<node> getComponent(index component);


	/**
	 * This method returns a map of component index
	 * to size of the component.
	 *
	 */
	std::map<index, count> getComponentSizes();

private:
	std::vector<node> component;
	std::vector<count> componentSizes;
};

}


#endif /* PARALLELCONNECTEDCOMPONENTS_H_ */
