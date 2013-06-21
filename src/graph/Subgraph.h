/*
 * Subgraph.h
 *
 *  Created on: Jun 13, 2013
 *      Author: forigem
 */

#ifndef SUBGRAPH_H_
#define SUBGRAPH_H_

#include <unordered_set>
#include <unordered_map>

#include "Graph.h"

namespace NetworKit {

class Subgraph {
public:
	Subgraph();
	virtual ~Subgraph();

	static Graph fromNodes(Graph G, std::unordered_set<node> nodeMap);

};


} /* namespace NetworKit */
#endif /* SUBGRAPH_H_ */
