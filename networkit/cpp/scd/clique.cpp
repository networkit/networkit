#include "clique.h"

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

	//TODO: use degeneracy ordering of vertices.
	for (const node& u : G.nodes()) {
		std::swap(pxvector[pxlookup[u]], pxvector[0]);
		std::swap(pxlookup[u], pxlookup[pxvector[0]]);

		uint32_t xcount = 0;
		uint32_t pcount = 0;
		G.forNeighborsOf(u, [&] (node v) {
			if (pxlookup[v] < xpbound) { // v is in X
				auto pxvec2 = pxvector[xpbound - xcount - 1];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound - xcount - 1]);
				std::swap(pxlookup[v], pxlookup[pxvec2]);
			} else { // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
				std::swap(pxlookup[v], pxlookup[pxvec2]);
			}
		});

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
	for (uint32_t i = xpbound; i < pbound; i++) {
		if (G.hasEdge(pxvector[i], u)) {
			break;
		}
		uint32_t xcount = 0, pcount = 0;
		G.forNeighborsOf(pxvector[i], [&] (node v) {
			if (pxlookup[v] < xpbound) { // v is in X
				auto pxvec2 = pxvector[xpbound - xcount - 1];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound - xcount - 1]);
				std::swap(pxlookup[v], pxlookup[pxvec2]);
			} else { // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
				std::swap(pxlookup[v], pxlookup[pxvec2]);
			}
		});
		std::vector<node> rplusv(r);
		rplusv.push_back(pxvector[i]);
		auto tmp = tomita(pxvector, pxlookup, neighbors, xpbound - xcount, xpbound, xpbound + pcount, rplusv);
		result.insert(result.end(), tmp.begin(), tmp.end());

		std::swap(pxvector[i], pxvector[xpbound]);
		std::swap(pxlookup[pxvector[i]], pxlookup[pxvector[xpbound]]);
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
	node maxnode;
	uint32_t maxval = 0;

	for (uint32_t i = xbound; i < pbound; i++) {
		auto inb = neighbors[i];
		uint32_t val = 0;
		/*do {
			val++;
		} while (pxlookup[inb[val]] >= xpbound && pxlookup[inb[val]] < pbound);
		*/
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

	return maxnode;
}

}
