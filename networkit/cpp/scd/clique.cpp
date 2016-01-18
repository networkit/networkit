#include "clique.h"

#include <cassert>
#include <algorithm>

namespace NetworKit {

Clique::Clique(const Graph& G) : G(G) {
}

//TODO: neighbors currently a dummy vector

std::vector<std::vector<node> > Clique::run() {
	std::vector<std::vector<node> > result;

	auto orderedNodes = getDegeneracyOrdering();
	
	std::vector<node> pxvector(G.numberOfNodes());
	std::unordered_map<node, uint32_t> pxlookup;
	std::vector<std::vector<node> > neighbors(pxvector.size());
	uint32_t ii = 0;
	for (const node u : orderedNodes) {
		pxvector[ii] = u;
		neighbors[ii] = G.neighbors(u);
		pxlookup[u] = ii;
		ii += 1;
	}

	#ifndef NDEBUG
	for (auto u : orderedNodes) {
		assert(std::find(pxvector.begin(), pxvector.end(), u) != pxvector.end());
	}
	#endif

	uint32_t xpbound = 1;
	for (const node& u : orderedNodes) {
		auto pxvec2 = pxvector[xpbound - 1];
		std::swap(pxvector[pxlookup[u]], pxvector[xpbound - 1]);
		//std::swap(pxlookup[u], pxlookup[pxvec2]);
		pxlookup[pxvec2] = pxlookup[u];
		pxlookup[u] = xpbound - 1;

		#ifndef NDEBUG
		for (auto v : orderedNodes) {
			if (v == u)
				break;

			assert(pxlookup[v] < xpbound);
		}
		#endif

		uint32_t xcount = 0;
		uint32_t pcount = 0;
		G.forNeighborsOf(u, [&] (node v) {

			#ifndef NDEBUG
			assert(pxlookup[v] >= 0);
			assert(pxlookup[v] < pxvector.size());

			assert(xcount <= xpbound);
			assert(pcount <= pxvector.size() - xpbound);
			#endif

			if (pxlookup[v] < xpbound) { // v is in X
				auto pxvec2 = pxvector[xpbound - xcount - 1];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound - xcount - 1]);
				//std::swap(pxlookup[v], pxlookup[pxvec2]);
				pxlookup[pxvec2] = pxlookup[v];
				pxlookup[v] = xpbound - xcount - 1;
				xcount += 1;
			} else { // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
				//std::swap(pxlookup[v], pxlookup[pxvec2]);
				pxlookup[pxvec2] = pxlookup[v];
				pxlookup[v] = xpbound + pcount;
				pcount += 1;
			}
		});

		#ifndef NDEBUG
		bool inRange = false;
		bool wasInRange = false;
		for (node v : pxvector) {
			if (G.hasEdge(u, v)) {
				assert(!wasInRange);
				inRange = true;
			} else {
				if (inRange) {
					wasInRange = true;
					inRange = false;
				}
			}
		}
		#endif

		std::vector<node> r;
		r.push_back(u);
		auto tmp = tomita(pxvector, pxlookup, neighbors, xpbound - xcount, xpbound, xpbound + pcount, r);
		result.insert(result.end(), tmp.begin(), tmp.end());
		xpbound += 1;
	}

	return result;
}

std::vector<std::vector<node> > Clique::tomita(std::vector<node>& pxvector, std::unordered_map<node, uint32_t>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound, std::vector<node>& r) {
	std::vector<std::vector<node> > result;
	if (xbound == pbound) { //if (X, P are empty)
		result.push_back(r);
		return result;
	}

	#ifndef NDEBUG
	assert(xbound >= 0);
	assert(xbound <= xpbound);
	assert(xpbound <= pbound);
	assert(pbound <= pxvector.size());
	#endif

	node u = findPivot(pxvector, pxlookup, neighbors, xbound, xpbound, pbound);
	std::vector<node> movedNodes;

	// this step is necessary as the next loop changes pxvector,
	// which prohibits iterating over it in the same loop.
	std::vector<node> toCheck;
	for (uint32_t i = xpbound; i < pbound; i++) {
		if (!G.hasEdge(pxvector[i], u)) {
			toCheck.push_back(pxvector[i]);
		}
	}

	for (auto pxveci : toCheck) {
		uint32_t xcount = 0, pcount = 0;
		G.forNeighborsOf(pxveci, [&] (node v) {
			if (pxlookup[v] < xpbound && pxlookup[v] >= xbound) { // v is in X
				auto pxvec2 = pxvector[xpbound - xcount - 1];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound - xcount - 1]);
				//std::swap(pxlookup[v], pxlookup[pxvec2]);
				pxlookup[pxvec2] = pxlookup[v];
				pxlookup[v] = xpbound - xcount - 1;
				xcount += 1;
			} else if (pxlookup[v] >= xpbound && pxlookup[v] < pbound){ // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
				//std::swap(pxlookup[v], pxlookup[pxvec2]);
				pxlookup[pxvec2] = pxlookup[v];
				pxlookup[v] = xpbound + pcount;
				pcount += 1;
			}
		});

		std::vector<node> rplusv(r);
		rplusv.push_back(pxveci);

		#ifndef NDEBUG
		assert(xpbound + pcount <= pbound);
		assert(xpbound - xcount >= xbound);
		#endif

		auto tmp = tomita(pxvector, pxlookup, neighbors, xpbound - xcount, xpbound, xpbound + pcount, rplusv);
		result.insert(result.end(), tmp.begin(), tmp.end());

		auto pxvec2 = pxvector[xpbound];
		std::swap(pxvector[pxlookup[pxveci]], pxvector[xpbound]);
		//std::swap(pxlookup[pxveci], pxlookup[pxvector[xpbound]]);
		pxlookup[pxvec2] = pxlookup[pxveci];
		pxlookup[pxveci] = xpbound;
		xpbound += 1;
		movedNodes.push_back(pxvector[xpbound - 1]);
	}
	
	for (node v : movedNodes) {
		//move from X -> P
		auto pxvec2 = pxvector[xpbound - 1];
		std::swap(pxvector[pxlookup[v]], pxvector[xpbound - 1]);
		//std::swap(pxlookup[v], pxlookup[pxvec2]);
		pxlookup[pxvec2] = pxlookup[v];
		pxlookup[v] = xpbound - 1;
		xpbound -= 1;
	}

	#ifndef NDEBUG
	for (node v : movedNodes) {
		assert(pxlookup[v] >= xpbound);
		assert(pxlookup[v] < pbound);
	}
	#endif

	return result;
}

node Clique::findPivot(std::vector<node>& pxvector, std::unordered_map<node, uint32_t>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound) {
	node maxnode = G.upperNodeIdBound() + 1;
	int32_t maxval = -1;

	for (uint32_t i = xbound; i < pbound; i++) {
		int32_t val = 0;
		G.forNeighborsOf(pxvector[i], [&] (node v) {
			if (pxlookup[v] >= xpbound && pxlookup[v] < pbound) {
				val++;
			}
		});

		if (val > maxval) {
			maxval = val;
			maxnode = pxvector[i];
		}
	}

	#ifndef NDEBUG
	assert(maxnode < G.upperNodeIdBound() + 1);
	#endif

	return maxnode;
}

std::vector<node> Clique::getDegeneracyOrdering() {
	std::vector<node> result;
	Graph dG(G);
	
	while (dG.numberOfNodes() > 0) {
		count minDegree = dG.numberOfNodes();
		node minNode = 0;
		dG.forNodes([&] (node u) {
			if (dG.degree(u) < minDegree) {
				minDegree = dG.degree(u);
				minNode = u;
			}
		});
		result.push_back(minNode);
		dG.forNeighborsOf(minNode, [&] (node v) {
			dG.removeEdge(minNode, v);
		});
		dG.removeNode(minNode);
	}

	return result;
}

}
