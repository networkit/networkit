/*
 * CommunityTrimming.h
 *
 * Created on: 10.06.2013
 * Author: cls, Yassine Marrakchi
 */

#ifndef COMMUNITYTRIMMING_H_
#define COMMUNITYTRIMMING_H_

#include <unordered_set>
#include "../graph/Graph.h"

namespace NetworKit {
/**
 * Trimming functions quantify how likely a node from the community shell
 * is to be removed from the community.
 */
class CommunityTrimming {
public:

	CommunityTrimming();

	virtual ~CommunityTrimming();

	virtual std::unordered_set<node> run(std::unordered_set<node>& community, const Graph& G) = 0;
};

/**
 * Remove nodes having more outgoing nodes than ingoing nodes
 */
class BoundarySharpness : public CommunityTrimming {
public:
	BoundarySharpness();

	virtual ~BoundarySharpness();

	virtual std::unordered_set<node> run(std::unordered_set<node>& community, const Graph& G);
};

/**
 * Return the same value for every node.
 *
 * This class is used for variants without trimming
 */
class DummyTrimming : public CommunityTrimming {
public:
	DummyTrimming();

	virtual ~DummyTrimming();

	virtual std::unordered_set<node> run(std::unordered_set<node>& community, const Graph& G);
};

} /* namespace NetworKit */
#endif /* COMMUNITYTRIMMING_H_ */
