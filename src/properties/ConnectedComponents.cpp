/*
 * ConnectedComponents.cpp
 *
 *  Created on: Sep 17, 2013
 *      Author: birderon
 */

#include "ConnectedComponents.h"

namespace NetworKit {

ConnectedComponents::ConnectedComponents() {

}

ConnectedComponents::ConnectedComponents(const Graph& g) {
	this->run(g);
}

ConnectedComponents::~ConnectedComponents() {

}

void ConnectedComponents::run(const Graph& g) {
	BFS bfs;
	count infDist = std::numeric_limits<count>::max();
	count n = g.numberOfNodes();
	this->assignedComponent = std::vector<node>(n,infDist);
	node startNode = 0;
	count idOfCC = 0;
	std::vector<count> resultOfBFS;
	//resultOfBFS = bfs.run(g,startNode);
	//for (count i = 0; i < resultOfBFS.size(); i++) {
	//	std::cout << resultOfBFS.at(i) << std::endl;
	//}
	do {
		//std::cout << "run: " << idOfCC << std::endl;
		resultOfBFS = bfs.run(g,startNode);
		startNode = infDist;
		count nodesInCC = 0;
		for (count i = 0; i < n;i++) {
			if (resultOfBFS.at(i) != infDist) {
				//TODO  what if the above and  'nodes[i] != infDist'?
				this->assignedComponent[i] = idOfCC;
				nodesInCC++;
			}
		}
		startNode = this->findStartNode();
		this->componentSizes.push_back(nodesInCC);
		idOfCC++;
	} while (startNode != infDist);

}

node ConnectedComponents::findStartNode() {
	count infDist = std::numeric_limits<count>::max();
	for (count i = 0; i < this->assignedComponent.size(); i++) {
		if (this->assignedComponent[i] == infDist) {
			return (node)i;
		}
	}
	return infDist;
}

std::vector<node> ConnectedComponents::getComponent(index component) {
	std::vector<node> nodesOfComponent;
	count i = 0;
	while (nodesOfComponent.size() < this->componentSizes.at(component)) {
		if (this->assignedComponent.at(i) == component) {
			nodesOfComponent.push_back(i);
		}
		i++;
	}
	return nodesOfComponent;
}

count ConnectedComponents::numberOfComponents() {
	return componentSizes.size();
}

count ConnectedComponents::sizeOfComponent(index component) {
	return componentSizes[component];
}

count ConnectedComponents::componentOfNode(node query) {
	return assignedComponent[query];
}

}

