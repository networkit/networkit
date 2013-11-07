/*
 * ConnectedComponents.h
 *
 *  Created on: Sep 17, 2013
 *      Author: Maximilian Vogel
 */

#ifndef CONNECTEDCOMPONENTS_H_
#define CONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../graph/BFS.h"

namespace NetworKit {

class ConnectedComponents {
public:
	ConnectedComponents();
	virtual ~ConnectedComponents();
	/**
	 * This method determines the connected components for the graph g.
	 *
	 * @param[in]	g	graph 
	 */
	void run(const Graph& g);

	/**
	 * This method returns the number of connected components.
	 */
	count numberOfComponents();

	/**
	 * This method returns the size of the given component.
	 *
	 * @param[in]	component	the index of questioned component
	 */
	count sizeOfComponent(index component);

	/**
	 * This method returns the the component in which node query is situated.
	 *
	 * @param[in]	query	the node whose component is asked for
	 */
	count componentOfNode(node query);

	/**
	 * This method returns a component consisting of a vector of nodes
	 *
	 * @param[in]	component	the index of the component that is asked for
	 */
	std::vector<node> getComponent(index component);

private:
	std::vector<node> assignedComponent;
	std::vector<count> componentSizes;
	node findStartNode();
};

}


#endif /* CONNECTEDCOMPONENTS_H_ */
