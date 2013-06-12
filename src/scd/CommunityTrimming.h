/*
 * CommunityTrimming.h
 *
 *  Created on: 10.06.2013
 *      Author: cls
 */

#ifndef COMMUNITYTRIMMING_H_
#define COMMUNITYTRIMMING_H_

#include <unordered_set>

#include "../graph/Graph.h"

namespace NetworKit {

class CommunityTrimming {
public:

	class NodeFitness {
		// TODO: constructor
		virtual double getValue(node v) = 0;
	};

	CommunityTrimming();

	virtual ~CommunityTrimming();

	virtual std::unordered_set<node> run(std::unordered_set<node>& community, Graph& G) = 0;
};

} /* namespace NetworKit */
#endif /* COMMUNITYTRIMMING_H_ */
