/*
 * GraphUpdater.h
 *
 *  Created on: 27.12.2013
 *      Author: cls
 */

#ifndef GRAPHUPDATER_H_
#define GRAPHUPDATER_H_

#include "../graph/Graph.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

class GraphUpdater {

public:

	GraphUpdater(Graph& G);

	void update(std::vector<GraphEvent>& stream);

private:

	Graph& G;
};

} /* namespace NetworKit */

#endif /* GRAPHUPDATER_H_ */
