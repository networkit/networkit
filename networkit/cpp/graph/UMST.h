#ifndef UMST_H
#define UMST_H

#include "Graph.h"
#include <limits>
#include "../structures/UnionFind.h"

namespace NetworKit {

/**
 * Union maximum spanning tree algorithm, computes the union of all maximum spanning trees.
 */
class UMST {
public:
	UMST(const Graph &G);
	
	Graph generate();

	template <typename A>
	inline Graph generate(const std::vector<A> &attribute);
	
	template <typename A>
	inline std::vector<bool> generateAttribute(const std::vector<A> &attribute);
private:
	template <typename A>
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

	template <typename A, typename F>
	inline void generate(std::vector< NetworKit::UMST::weightedEdge< A > > weightedEdges, F callback);
	
	template <typename A>
	inline std::vector<weightedEdge<A> > genWeightedEdges(const std::vector<A> &attribute);

	const Graph &G;
};

template <typename A>
std::vector< UMST::weightedEdge< A > > UMST::genWeightedEdges(const std::vector< A > &attribute) {
	std::vector<weightedEdge<A> > weightedEdges;

	if (!G.hasEdgeIds()) {
		throw std::runtime_error("Error: Edges of G must be indexed for using edge attributes");
	}

	weightedEdges.reserve(G.numberOfEdges());

	G.forEdges([&](node u, node v, edgeid eid) {
		weightedEdges.emplace_back(u, v, attribute[eid], eid);
	});

	return weightedEdges;
}


template <typename A>
Graph UMST::generate(const std::vector< A > &attribute) {
	Graph result(G.copyNodes());

	generate(genWeightedEdges(attribute), [&](weightedEdge<A> &e) {
		result.addEdge(e.u, e.v);
	});

	return result;
}

template <typename A>
std::vector< bool > UMST::generateAttribute(const std::vector< A > &attribute) {
	std::vector<bool> result(G.upperEdgeIdBound(), false);

	generate(genWeightedEdges(attribute), [&](weightedEdge<A> &e) {
		result[e.eid] = true;
	});

	return result;
}


template <typename A, typename F>
void UMST::generate(std::vector< UMST::weightedEdge< A > > weightedEdges, F callback) {
	std::sort(weightedEdges.begin(), weightedEdges.end(), std::greater<weightedEdge<A> >());

	A currentAttribute = std::numeric_limits<A>::max();

	std::vector<std::pair<node, node> > nodesToMerge;
	UnionFind uf(G.upperNodeIdBound());

	for (weightedEdge<A> e : weightedEdges) {
		if (e.attribute != currentAttribute) {
			for (auto candidate : nodesToMerge) {
				uf.merge(candidate.first, candidate.second);
			}

			nodesToMerge.clear();
			currentAttribute = e.attribute;
		}

		
		if (uf.find(e.u) != uf.find(e.v)) {
			callback(e);
			nodesToMerge.emplace_back(e.u, e.v);
	
		}
	}
}


} // namespace NetworKit

#endif // UMST_H
