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

	void update(Graph& G, std::vector<GraphEvent>& stream);
};

} /* namespace NetworKit */

#endif /* GRAPHUPDATER_H_ */
