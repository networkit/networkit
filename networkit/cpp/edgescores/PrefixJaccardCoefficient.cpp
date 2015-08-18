/*
 *
 */

#include "PrefixJaccardCoefficient.h"


namespace NetworKit {

template <typename AttributeT>
PrefixJaccardCoefficient<AttributeT>::PrefixJaccardCoefficient(const Graph &G, const std::vector< AttributeT > &attribute) : G(G), inAttribute(attribute) {
}

template <typename AttributeT>
void PrefixJaccardCoefficient<AttributeT>::run() {
	if (!G.hasEdgeIds()) throw std::runtime_error("Error, edges need to be indexed first");

	outAttribute.clear();
	outAttribute.resize(G.upperEdgeIdBound());

	struct RankedEdge {
		node u;
		AttributeT att;
		count rank;

		RankedEdge(node u, AttributeT att, count rank) : u(u), att(att), rank(rank) {};

		bool operator<(const RankedEdge &other) const {
			return std::tie(rank, att, u) < std::tie(other.rank, other.att, other.u);
		};

		bool operator>(const RankedEdge &other) const {
			return std::tie(rank, att, u) > std::tie(other.rank, other.att, other.u);
		};
	};

	std::vector<std::vector<RankedEdge> > rankedEdges(G.upperNodeIdBound());

	G.balancedParallelForNodes([&](node u) {
		if (G.degree(u) == 0) return;

		rankedEdges[u].reserve(G.degree(u));

		G.forEdgesOf(u, [&](node _u, node w, edgeid eid) {
			rankedEdges[u].emplace_back(w, inAttribute[eid], 0);
		});

		std::sort(rankedEdges[u].begin(), rankedEdges[u].end(), std::greater<RankedEdge>());

		AttributeT curVal = rankedEdges[u].front().att;
		count curRank = 0;
		count numEqual = 0;
		for (RankedEdge& cur : rankedEdges[u]) {
			if (curVal != cur.att) {
				curRank += numEqual;
				curVal = cur.att;
				numEqual = 1;
			} else {
				++numEqual;
			}

			cur.rank = curRank;
		}
	});

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		std::set<node> uNeighbors, vNeighbors;
		count curRank = 0;
		double bestJaccard = 0;

		auto uIt = rankedEdges[u].begin();
		auto vIt = rankedEdges[v].begin();
		count commonNeighbors = 0;

		while (uIt != rankedEdges[u].end() || vIt != rankedEdges[v].end()) {
			while (uIt != rankedEdges[u].end() && curRank == uIt->rank) {
				if (uIt->u == v) {
					++uIt;
					continue;
				}

				if (vNeighbors.erase(uIt->u)) {
					++commonNeighbors;
				} else {
					uNeighbors.insert(uIt->u);
				}

				++uIt;
			}

			while (vIt != rankedEdges[v].end() && curRank == vIt->rank) {
				if (vIt->u == u) {
					++vIt;
					continue;
				}

				if (uNeighbors.erase(vIt->u)) {
					++commonNeighbors;
				} else {
					vNeighbors.insert(vIt->u);
				}

				++vIt;
			}

			bestJaccard = std::max(bestJaccard, commonNeighbors * 1.0 / (uNeighbors.size() + vNeighbors.size() + commonNeighbors));

			++curRank;
		}

		outAttribute[eid] = bestJaccard;
	});

	hasRun = true;
}


template <typename AttributeT>
std::vector< double > PrefixJaccardCoefficient<AttributeT>::getAttribute() {
	if (!hasRun) throw std::runtime_error("Error, run must be called first");

	return outAttribute;
}


template class PrefixJaccardCoefficient<double>;
template class PrefixJaccardCoefficient<count>;

}