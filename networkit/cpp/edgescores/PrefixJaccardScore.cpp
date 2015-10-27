/*
 * PrefixJaccardScore.h
 *
 *  Created on: 26.09.2014
 *      Author: Michael Hamann
 */

#include "PrefixJaccardScore.h"

namespace NetworKit {

template <typename AttributeT>
PrefixJaccardScore<AttributeT>::PrefixJaccardScore(const Graph &G, const std::vector< AttributeT > &attribute) :
	EdgeScore<double>(G), inAttribute(attribute) {
}

template <typename AttributeT>
void PrefixJaccardScore<AttributeT>::run() {
	//this-> required to access members of base class, since this is a template class.
	if (!this->G.hasEdgeIds()) throw std::runtime_error("Error, edges need to be indexed first");

	this->scoreData.clear();
	this->scoreData.resize(this->G.upperEdgeIdBound());

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

	std::vector<std::vector<RankedEdge> > rankedEdges(this->G.upperNodeIdBound());

	this->G.balancedParallelForNodes([&](node u) {
		if (this->G.degree(u) == 0) return;

		rankedEdges[u].reserve(this->G.degree(u));

		this->G.forEdgesOf(u, [&](node _u, node w, edgeid eid) {
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

	this->G.parallelForEdges([&](node u, node v, edgeid eid) {
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

		this->scoreData[eid] = bestJaccard;
	});

	this->hasRun = true;
}

template <typename AttributeT>
double PrefixJaccardScore<AttributeT>::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

template <typename AttributeT>
double PrefixJaccardScore<AttributeT>::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

template class PrefixJaccardScore<double>;
template class PrefixJaccardScore<count>;

}
