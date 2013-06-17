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

	CommunityTrimming();

	virtual ~CommunityTrimming();

	virtual std::unordered_set<node> run(std::unordered_set<node>& community, Graph& G) = 0;
};

class BoundarySharpness : public CommunityTrimming {
public:
	BoundarySharpness();

	virtual ~BoundarySharpness();

	virtual std::unordered_set<node> run(std::unordered_set<node>& community, Graph& G);
};

class DummyTrimming : public CommunityTrimming {
public:
	DummyTrimming();

	virtual ~DummyTrimming();

	virtual std::unordered_set<node> run(std::unordered_set<node>& community, Graph& G);
};

} /* namespace NetworKit */
#endif /* COMMUNITYTRIMMING_H_ */
