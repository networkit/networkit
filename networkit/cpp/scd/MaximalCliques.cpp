#include "MaximalCliques.h"

#include <cassert>
#include <algorithm>

namespace {
	using NetworKit::node;
	using NetworKit::index;

	struct SwapFunctor {
		std::vector<node> &pxvector;
		std::vector<node> &pxlookup;

		SwapFunctor(std::vector<node> &pxvector, std::vector<index>& pxlookup) : pxvector(pxvector), pxlookup(pxlookup) {
		}

		void operator()(node u, index pos) {
			node pxvec2 = pxvector[pos];
			std::swap(pxvector[pxlookup[u]], pxvector[pos]);
			pxlookup[pxvec2] = pxlookup[u];
			pxlookup[u] = pos;
		}
	};
}

namespace NetworKit {

MaximalCliques::MaximalCliques(const Graph& G) : G(G) {
}

//TODO: neighbors currently a dummy vector

std::vector<std::vector<node> > MaximalCliques::run() {
	std::vector<std::vector<node> > result;

	auto orderedNodes = getDegeneracyOrdering();

	std::vector<node> pxvector(G.numberOfNodes());
	std::vector<index> pxlookup(G.upperNodeIdBound());
	
	SwapFunctor swapNodeToPos(pxvector, pxlookup);

	uint32_t ii = 0;
	for (const node u : orderedNodes) {
		pxvector[ii] = u;
		pxlookup[u] = ii;
		ii += 1;
	}

	#ifndef NDEBUG
	for (auto u : orderedNodes) {
		assert(pxvector[pxlookup[u]] == u);
	}
	#endif

	// Store out-going neighbors in the direction of higher core numbers.
	// This means that the out-degree is bounded by the maximum core number.
	const StaticOutGraph outGraph(G, pxlookup);

	uint32_t xpbound = 1;
	for (const node& u : orderedNodes) {
		swapNodeToPos(u, xpbound-1);

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
				swapNodeToPos(v, xpbound - xcount - 1);
				xcount += 1;
			} else { // v is in P
				swapNodeToPos(v, xpbound + pcount);
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

		std::vector<node> r = {u};
		tomita(outGraph, pxvector, pxlookup, xpbound - xcount, xpbound, xpbound + pcount, r, result);

		xpbound += 1;
	}

	return result;
}

void MaximalCliques::tomita(const StaticOutGraph& outGraph, std::vector<node>& pxvector, std::vector<index>& pxlookup, uint32_t xbound, uint32_t xpbound, uint32_t pbound, std::vector<node>& r, std::vector<std::vector<node>>& result) {
	if (xbound == pbound) { //if (X, P are empty)
		result.push_back(r);
		return;
	}

	if (xpbound == pbound) return;

	SwapFunctor swapNodeToPos(pxvector, pxlookup);

	#ifndef NDEBUG
	assert(xbound >= 0);
	assert(xbound <= xpbound);
	assert(xpbound <= pbound);
	assert(pbound <= pxvector.size());
	#endif

	node u = findPivot(outGraph, pxvector, pxlookup, xbound, xpbound, pbound);
	std::vector<node> movedNodes;

	// Find all nodes in P that are not neighbors of the pivot
	// this step is necessary as the next loop changes pxvector,
	// which prohibits iterating over it in the same loop.
	std::vector<node> toCheck;

	// Step 1: mark all outgoing neighbors of the pivot in P
	std::vector<bool> pivotNeighbors(pbound - xpbound);
	outGraph.forOutEdgesOf(u, [&](node v) {
		index vpos = pxlookup[v];
		if (vpos >= xpbound && vpos < pbound) {
			pivotNeighbors[vpos - xpbound] = true;
		}
	});

	// Step 2: for all not-yet marked notes check if they have the pivot as neighbor.
	// If not: they are definitely a non-neighbor.
	for (index i = xpbound; i < pbound; i++) {
		if (!pivotNeighbors[i - xpbound]) {
			node p = pxvector[i];

			if (!outGraph.hasNeighbor(p, u)) {
				toCheck.push_back(p);
			}
		}
	}

	for (auto pxveci : toCheck) {
		uint32_t xcount = 0, pcount = 0;

		// Group all neighbors of pxveci in P \cup X around xpbound.
		// Step 1: collect all outgoing neighbors of pxveci
		outGraph.forOutEdgesOf(pxveci, [&](node v) {
			if (pxlookup[v] < xpbound && pxlookup[v] >= xbound) { // v is in X
				swapNodeToPos(v, xpbound - xcount - 1);
				xcount += 1;
			} else if (pxlookup[v] >= xpbound && pxlookup[v] < pbound){ // v is in P
				swapNodeToPos(v, xpbound + pcount);
				pcount += 1;
			}
		});

		// Step 2: collect all nodes in X that have not yet been collected
		// and that have pxveci as outgoing neighbor.
		for (index i = xbound; i < xpbound;) {
			// stop if we have reached the collected neighbors
			if (i == xpbound - xcount) break;
			node x = pxvector[i];

			if (outGraph.hasNeighbor(x, pxveci)) {
				swapNodeToPos(x, xpbound - xcount - 1);
				xcount += 1;
			} else {
				// Advance only if we did not swap otherwise we have already
				// a next candidate at position i.
				++i;
			}
		}

		// Step 3: collect all nodes in P that have not yet been collected
		// and that have pxveci as outgoing neighbor.
		for (index i = xpbound + pcount; i < pbound; ++i) {
			node p = pxvector[i];

			if (outGraph.hasNeighbor(p, pxveci)) {
				swapNodeToPos(p, xpbound + pcount);
				pcount += 1;
			}
		}

		r.push_back(pxveci);

		#ifndef NDEBUG
		assert(xpbound + pcount <= pbound);
		assert(xpbound - xcount >= xbound);
		#endif

		tomita(outGraph, pxvector, pxlookup, xpbound - xcount, xpbound, xpbound + pcount, r, result);

		r.pop_back();

		swapNodeToPos(pxveci, xpbound);
		xpbound += 1;
		assert(pxvector[xpbound - 1] == pxveci);
		movedNodes.push_back(pxveci);
	}

	for (node v : movedNodes) {
		//move from X -> P
		swapNodeToPos(v, xpbound - 1);
		xpbound -= 1;
	}

	#ifndef NDEBUG
	for (node v : movedNodes) {
		assert(pxlookup[v] >= xpbound);
		assert(pxlookup[v] < pbound);
	}
	#endif
}

node MaximalCliques::findPivot(const StaticOutGraph& outGraph, std::vector<node>& pxvector, std::vector<index>& pxlookup, uint32_t xbound, uint32_t xpbound, uint32_t pbound) {
	node maxnode = none;
	int32_t maxval = -1;

	// Counts for every node in X \cup P how many outgoing neighbors it has in P
	std::vector<count> pivotNeighbors(pbound - xbound);
	const count psize = pbound-xpbound;

	// Step 1: for all nodes in X count how many outgoing neighbors they have in P
	for (index i = 0; i < xpbound - xbound; i++) {
		node u = pxvector[i + xbound];
		outGraph.forOutEdgesOf(u, [&](node v) {
			if (pxlookup[v] >= xpbound && pxlookup[v] < pbound) {
				++pivotNeighbors[i];
			}
		});

		// If a node has |P| neighbors, we cannot find a better candidate
		if (pivotNeighbors[i] == psize) return u;
	}

	// Step 2: for all nodes in P
	// a) increase counts for every neighbor in P \cup X to account for incoming neighbors
	// b) count all outgoing neighbors in P
	for (index i = xpbound - xbound; i < pivotNeighbors.size(); ++i) {
		node u = pxvector[i + xbound];
		outGraph.forOutEdgesOf(u, [&](node v) {
			index neighborPos = pxlookup[v];
			if (neighborPos >= xbound && neighborPos < pbound) {
				++pivotNeighbors[neighborPos-xbound];

				if (neighborPos >= xpbound) {
					++pivotNeighbors[i];
				}
			}
		});
	}

	// Step 3: find maximum
	for (index i = 0; i < pivotNeighbors.size(); ++i) {
		if (static_cast<int64_t>(pivotNeighbors[i]) > maxval) {
			maxval = pivotNeighbors[i];
			maxnode = pxvector[i + xbound];
		}
	}

	#ifndef NDEBUG
	assert(maxnode < G.upperNodeIdBound() + 1);
	#endif

	return maxnode;
}

std::vector<node> MaximalCliques::getDegeneracyOrdering() {
	std::vector<node> result;
	result.reserve(G.numberOfNodes());

	/* Main data structure: buckets of nodes indexed by their remaining degree. */
	index z = G.upperNodeIdBound();
	std::vector<node> queue(G.numberOfNodes());
	std::vector<index> nodePtr(z);
	std::vector<index> degreeBegin(G.numberOfNodes());
	std::vector<count> degree(z);       // tracks degree during algo

	/* Bucket sort  by degree */
	/* 1) bucket sizes */
	G.forNodes([&](node u) {
		count deg = G.degree(u);
		degree[u] = deg;
		++degreeBegin[deg];
	});

	index sum = 0; // 2) exclusive in-place prefix sum
	for (index i = 0; i < degreeBegin.size(); ++i) {
		index tmp = degreeBegin[i];
		degreeBegin[i] = sum;
		sum += tmp;
	}

	/* 3) Sort nodes/place them in queue */
	G.forNodes([&](node u) {
		count deg = degree[u];
		index pos = degreeBegin[deg];
		++degreeBegin[deg];
		queue[pos] = u;
		nodePtr[u] = pos;
	});

	/* 4) restore exclusive prefix sum */
	index tmp = 0; // move all one forward
	for (index i = 0; i < degreeBegin.size(); ++i) {
		std::swap(tmp, degreeBegin[i]);
	}

	/* Current core and and computed scoreData values. */
	index core = 0;

	/* Main loop: Successively "remove" nodes by setting them not alive after processing them. */
	for (index i = 0; i < G.numberOfNodes(); ++i) {
		node u = queue[i];
		core = std::max(core, degree[u]); // core is maximum of all previously seen degrees

		result.emplace_back(u);

		/* Remove a neighbor by decreasing its degree and changing its position in the queue */
		auto removeNeighbor = [&](node v) {
			if (nodePtr[v] > i) { // only nodes that are after the current node need to be considered
				// adjust the degree
				count oldDeg = degree[v];
				--degree[v];
				count newDeg = oldDeg - 1;

				// Degrees smaller than the current degree can be before the current position
				// Correct those that we need. Note that as we decrease degrees only by one
				// all degreeBegin values larger than oldDeg will have been adjusted already.
				if (degreeBegin[oldDeg] <= i) {
					degreeBegin[oldDeg] = i + 1;
				}

				if (degreeBegin[newDeg] <= i) {
					degreeBegin[newDeg] = i + 1;
				}

				// Swap v with beginning of the current bucket.
				index oldPos = nodePtr[v];
				index newPos = degreeBegin[oldDeg];
				node nodeToSwap = queue[newPos];
				std::swap(queue[oldPos], queue[newPos]);
				std::swap(nodePtr[nodeToSwap], nodePtr[v]);

				// Move bucket border, v is now in the previous bucket, i.e. the bucket of its new degree
				++degreeBegin[oldDeg];
			}
		};

		/* Remove u and its incident edges. */
		G.forNeighborsOf(u, removeNeighbor);
	}
	
	return result;
}

}
