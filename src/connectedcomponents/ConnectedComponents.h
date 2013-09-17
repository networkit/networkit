/*
 * ConnectedComponents.h
 *
 *  Created on: Sep 17, 2013
 *      Author: birderon
 */

#ifndef CONNECTEDCOMPONENTS_H_
#define CONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../graph/BFS.h"

namespace NetworKit {

class ConnectedComponents {
public:
	ConnectedComponents();
	ConnectedComponents(const Graph& g);
	virtual ~ConnectedComponents();
	void run(const Graph& g);
	count numberOfComponents();
	count sizeOfComponent(index component);
	count componentOfNode(node query);
	std::vector<node> getComponent(index component);

private:
	std::vector<node> assignedComponent;
	std::vector<count> componentSizes;
	node findStartNode();
};

}


#endif /* CONNECTEDCOMPONENTS_H_ */
