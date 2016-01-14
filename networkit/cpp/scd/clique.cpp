#include "clique.h"

#include <iostream>
#include <cassert>

namespace NetworKit {

Clique::Clique(const Graph& G) : G(G) {
}

//TODO: neighbors currently a dummy vector

std::vector<std::vector<node> > Clique::run(node& seed) {
	std::vector<std::vector<node> > result;

	std::vector<node> pxvector = G.nodes();
	//if there are large gaps between the bounds, change to unordered_map
	std::vector<node> pxlookup(G.upperNodeIdBound());
	G.forNodes([&] (node u) { 
		pxlookup[u] = std::find(pxvector.begin(), pxvector.end(), u) - pxvector.begin();
	});
	std::vector<std::vector<node> > neighbors(pxvector.size());
	G.forNodes([&] (node u) {
		neighbors[pxlookup[u]] = G.neighbors(u);
	});

	uint32_t xpbound = 1;

	auto orderedNodes = getDegeneracyOrdering(G.nodes());
	for (const node& u : orderedNodes) {
		auto pxvec2 = pxvector[xpbound - 1];
		std::swap(pxvector[pxlookup[u]], pxvector[xpbound - 1]);
		std::swap(pxlookup[u], pxlookup[pxvec2]);

		uint32_t xcount = 0;
		uint32_t pcount = 0;
		G.forNeighborsOf(u, [&] (node v) {
			if (pxlookup[v] < xpbound) { // v is in X
				auto pxvec2 = pxvector[xpbound - xcount - 1];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound - xcount - 1]);
				std::swap(pxlookup[v], pxlookup[pxvec2]);
				xcount += 1;
			} else { // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
				std::swap(pxlookup[v], pxlookup[pxvec2]);
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
		xpbound++;
	}

	return result;
}

std::vector<std::vector<node> > Clique::tomita(std::vector<node>& pxvector, std::vector<node>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound, std::vector<node>& r) {
	std::vector<std::vector<node> > result;
	if (xbound == pbound) { //if (X, P are empty)
		result.push_back(r);
		return result;
	}

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
				std::swap(pxlookup[v], pxlookup[pxvec2]);
				xcount += 1;
			} else if (pxlookup[v] >= xpbound && pxlookup[v] < pbound){ // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
				std::swap(pxlookup[v], pxlookup[pxvec2]);
				pcount += 1;
			}
		});

		// add pxvector[i] to r & move out P
		std::vector<node> rplusv(r);
		rplusv.push_back(pxveci);

		#ifndef NDEBUG
		assert(xpbound + pcount <= pbound);
		#endif

		auto tmp = tomita(pxvector, pxlookup, neighbors, xpbound - xcount, xpbound, xpbound + pcount, rplusv);
		result.insert(result.end(), tmp.begin(), tmp.end());

		std::swap(pxvector[pxlookup[pxveci]], pxvector[xpbound]);
		std::swap(pxlookup[pxveci], pxlookup[pxvector[xpbound]]);
		xpbound += 1;
		movedNodes.push_back(pxvector[xpbound - 1]);
	}
	
	for (node v : movedNodes) {
		//move from X -> P
		auto pxvec2 = pxvector[xpbound - 1];
		std::swap(pxvector[pxlookup[v]], pxvector[xpbound - 1]);
		std::swap(pxlookup[v], pxlookup[pxvec2]);
	}

	return result;
}

node Clique::findPivot(std::vector<node>& pxvector, std::vector<node>& pxlookup, std::vector<std::vector<node> >& neighbors, uint32_t xbound, uint32_t xpbound, uint32_t pbound) {
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

std::vector<node> Clique::getDegeneracyOrdering(std::vector<node> nodes) {
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
