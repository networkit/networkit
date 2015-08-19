#ifndef MST_H
#define MST_H

#include "Graph.h"
#include <limits>
#include "../structures/UnionFind.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

/**
 * Union maximum spanning tree algorithm, computes the union of all maximum spanning trees.
 */
class MST {
public:
	MST(const Graph &G);

	MST(const Graph &G, const std::vector<double> &attribute);

	void run();

	std::vector<bool> getAttribute(bool move = false);
	bool inMST(node u, node v) const;
	bool inMST(edgeid eid) const;

	Graph getMST(bool move = false);

private:
	struct weightedEdge {
		double attribute;
		node u;
		node v;
		edgeid eid;
		index rand;

		bool operator>(const weightedEdge &other) const {
			return
			(attribute > other.attribute)
			||
			(attribute == other.attribute &&
				(rand > other.rand ||
				(rand == other.rand && (u > other.u || (u == other.u && v > other.v)))));
		};
		weightedEdge(node u, node v, double attribute, edgeid eid = 0) : attribute(attribute), u(u), v(v), eid(eid), rand(Aux::Random::integer()) {};
	};

	const Graph &G;
	std::vector<weightedEdge> weightedEdges;

	Graph mst;
	std::vector<bool> mstAttribute;

	bool hasWeightedEdges;
	bool hasMST;
	bool hasAttribute;
};



} // namespace NetworKit

#endif // MST_H
