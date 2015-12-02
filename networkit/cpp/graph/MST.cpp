/*
 *
 */

#include "MST.h"

namespace NetworKit {

MST::MST(const Graph &G) : G(G), hasWeightedEdges(false), hasMST(false), hasAttribute(false) { };

MST::MST(const Graph &G, const std::vector< double > &attribute) : G(G), hasWeightedEdges(false), hasMST(false), hasAttribute(false) {
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

void MST::run() {
	mst = G.copyNodes();

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
		mstAttribute.resize(G.upperEdgeIdBound(), false);
		hasAttribute = true;
	}

	std::sort(weightedEdges.begin(), weightedEdges.end(), std::greater<weightedEdge>());

	UnionFind uf(G.upperNodeIdBound());

	for (weightedEdge e : weightedEdges) {
		if (uf.find(e.u) != uf.find(e.v)) {
			if (useEdgeWeights) {
				mst.addEdge(e.u, e.v, e.attribute);
			} else {
				mst.addEdge(e.u, e.v);
			}

			if (hasAttribute) {
				mstAttribute[e.eid] = true;
			}

			uf.merge(e.u, e.v);
		}
	}

	hasMST = true;
}

bool MST::inMST(edgeid eid) const {
	if (!hasAttribute) throw std::runtime_error("Error: Either the attribute hasn't be calculated yet or the graph has no edge ids.");

	return mstAttribute[eid];
}

bool MST::inMST(node u, node v) const {
	if (hasMST) {
		return mst.hasEdge(u, v);
	} else if (hasAttribute) {
		return mstAttribute[G.edgeId(u, v)];
	} else {
		throw std::runtime_error("Error: The run() method must be executed first");
	}
}

std::vector< bool > MST::getAttribute(bool move) {
	std::vector<bool> result;

	if (!hasAttribute) throw std::runtime_error("Error: The run() method must be executed first");

	if (move) {
		result = std::move(mstAttribute);
		hasAttribute = false;
	} else {
		result = mstAttribute;
	}

	return result;
}

Graph MST::getMST(bool move) {
	Graph result;

	if (!hasMST) throw std::runtime_error("Error: The run() method must be executed first");

	if (move) {
		result = std::move(mst);
		hasMST = false;
	} else {
		result = mst;
	}

	return result;
}

} // namespace NetworKit