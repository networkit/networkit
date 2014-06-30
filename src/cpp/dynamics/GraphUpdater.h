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

/**
 * @ingroup dynamics
 */
class GraphUpdater {

public:

	GraphUpdater(Graph& G);

	void update(std::vector<GraphEvent>& stream);

	std::vector<std::pair<count, count> > getSizeTimeline();

private:

	Graph& G;
	std::vector<std::pair<count, count> > size;
};

} /* namespace NetworKit */

#endif /* GRAPHUPDATER_H_ */
