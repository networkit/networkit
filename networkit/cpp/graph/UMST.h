#ifndef UMST_H
#define UMST_H

#include "Graph.h"
#include <limits>
#include "../structures/UnionFind.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

/**
 * Union maximum spanning tree algorithm, computes the union of all maximum spanning trees.
 */
template <typename A = edgeweight>
class UMST {
public:
	UMST(const Graph &G);

	UMST(const Graph &G, const std::vector<A> &attribute);
	
	void run();

	std::vector<bool> getAttribute(bool move = false);
	bool inUMST(node u, node v) const;
	bool inUMST(edgeid eid) const;

	Graph getUMST(bool move = false);

private:
	struct weightedEdge {
		A attribute;
		node u;
		node v;
		edgeid eid;
		
		bool operator>(const weightedEdge &other) const {
			return (attribute > other.attribute) || (attribute == other.attribute && (u > other.u || (u == other.u && v > other.v)));
		};
		weightedEdge(node u, node v, A attribute, edgeid eid = 0) : attribute(attribute), u(u), v(v), eid(eid) {};
	};

	const Graph &G;
	std::vector<weightedEdge> weightedEdges;

	Graph umst;
	std::vector<bool> umstAttribute;

	bool hasWeightedEdges;
	bool hasUMST;
	bool hasAttribute;
};

template <typename A>
UMST<A>::UMST(const Graph &G) : G(G), hasWeightedEdges(false), hasUMST(false), hasAttribute(false) { };

template <typename A>
UMST<A>::UMST(const Graph &G, const std::vector< A > &attribute) : G(G), hasWeightedEdges(false), hasUMST(false), hasAttribute(false) {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("Error: Edges of G must be indexed for using edge attributes");
	}

	weightedEdges.reserve(G.numberOfEdges());

	G.forEdges([&](node u, node v, edgeid eid) {
		weightedEdges.emplace_back(u, v, attribute[eid], eid);
	});

	INFO(weightedEdges.size(), " weighted edges saved");

	hasWeightedEdges = true;
}

template <typename A>
void UMST<A>::run() {
	umst = G.copyNodes();

	bool useEdgeWeights = false;

	if (!hasWeightedEdges) {
		weightedEdges.reserve(G.numberOfEdges());

		G.forEdges([&](node u, node v, edgeweight weight, edgeid eid) {
			weightedEdges.emplace_back(u, v, weight, eid);
		});

		hasWeightedEdges = true;
		useEdgeWeights = true;
	}

	if (G.hasEdgeIds()) {
		umstAttribute.resize(G.upperEdgeIdBound(), false);
		hasAttribute = true;
	}

	std::sort(weightedEdges.begin(), weightedEdges.end(), std::greater<weightedEdge>());

	A currentAttribute = std::numeric_limits<A>::max();

	std::vector<std::pair<node, node> > nodesToMerge;
	UnionFind uf(G.upperNodeIdBound());

	for (weightedEdge e : weightedEdges) {
		if (e.attribute != currentAttribute) {
			for (auto candidate : nodesToMerge) {
				uf.merge(candidate.first, candidate.second);
			}

			nodesToMerge.clear();
			currentAttribute = e.attribute;
		}


		if (uf.find(e.u) != uf.find(e.v)) {
			if (useEdgeWeights) {
				umst.addEdge(e.u, e.v, e.attribute);
			} else {
				umst.addEdge(e.u, e.v);
			}

			if (hasAttribute) {
				umstAttribute[e.eid] = true;
			}

			nodesToMerge.emplace_back(e.u, e.v);

		}
	}

	hasUMST = true;
}

template <typename A>
bool UMST<A>::inUMST(edgeid eid) const {
	if (!hasAttribute) throw std::runtime_error("Error: Either the attribute hasn't be calculated yet or the graph has no edge ids.");

	return umstAttribute[eid];
}

template <typename A>
bool UMST<A>::inUMST(node u, node v) const {
	if (hasUMST) {
		return umst.hasEdge(u, v);
	} else if (hasAttribute) {
		return umstAttribute[G.edgeId(u, v)];
	} else {
		throw std::runtime_error("Error: The run() method must be executed first");
	}
}

template <typename A>
std::vector< bool > UMST<A>::getAttribute(bool move) {
	std::vector<bool> result;

	if (!hasAttribute) throw std::runtime_error("Error: The run() method must be executed first");

	if (move) {
		result = std::move(umstAttribute);
		hasAttribute = false;
	} else {
		result = umstAttribute;
	}

	return result;
}

template <typename A>
Graph UMST<A>::getUMST(bool move) {
	Graph result;

	if (!hasUMST) throw std::runtime_error("Error: The run() method must be executed first");

	if (move) {
		result = std::move(umst);
		hasUMST = false;
	} else {
		result = umst;
	}

	return result;
}


} // namespace NetworKit

#endif // UMST_H
