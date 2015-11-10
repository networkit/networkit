/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#ifndef CONNECTEDCOMPONENTS_H_
#define CONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../graph/BFS.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"
#include <unordered_set>

namespace NetworKit {

/**
 * @ingroup components
 * Determines the connected components of an undirected graph.
 */
class ConnectedComponents : public Algorithm {
public:
	/**
	 * Create ConnectedComponents class for Graph @a G.
	 *
	 * @param G The graph.
	 */
	ConnectedComponents(const Graph& G);

	/**
	 * This method determines the connected components for the graph given in the constructor.
	 */
	void run();

	/**
	 * Get the number of connected components.
	 *
	 * @return The number of connected components.
	 */
	count numberOfComponents();

	/**
	 * Get the the component in which node @a u is situated.
	 *
	 * @param[in]	u	The node whose component is asked for.
	 */
	count componentOfNode(node u);


	/**
	 * Get a Partition that represents the components.
	 *
	 * @return A partition representing the found components.
	 */
	Partition getPartition();

    /**
     *Return the map from component to size
     */
    std::map<index, count> getComponentSizes();

    /**
     * @return Vector of components, each stored as (unordered) set of nodes.
     */
    std::vector<std::vector<node> > getComponents();


private:
	const Graph& G;
	Partition component;
	count numComponents;
	bool hasRun;
};

inline count ConnectedComponents::componentOfNode(node u) {
	assert (component[u] != none);
	if (!hasRun) throw std::runtime_error("run method has not been called");
	return component[u];
}

inline count ConnectedComponents::numberOfComponents() {
	if (!hasRun) throw std::runtime_error("run method has not been called");
	return this->numComponents;
}

}


#endif /* CONNECTEDCOMPONENTS_H_ */
