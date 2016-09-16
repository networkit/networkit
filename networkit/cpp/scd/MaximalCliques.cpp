#include "MaximalCliques.h"

#include <cassert>
#include <algorithm>

namespace {
#include <string>
#include <chrono>
#include <iostream>

class ScopedTimer { 
	private:
		std::chrono::time_point<std::chrono::high_resolution_clock> start;
		std::string message;
	public:
		ScopedTimer(const std::string& message) : start(std::chrono::high_resolution_clock::now()), message(message) {};
		~ScopedTimer() {
			std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now()-start;
			INFO(message, " needed: ", elapsed_seconds.count(), "s");

		};
};
};

namespace NetworKit {

MaximalCliques::MaximalCliques(const Graph& G) : G(G) {
}

//TODO: neighbors currently a dummy vector

std::vector<std::vector<node> > MaximalCliques::run() {
	std::vector<std::vector<node> > result;

	ScopedTimer *timer = new ScopedTimer("Degeneracy ordering");

	auto orderedNodes = getDegeneracyOrdering();

	delete timer;
	timer = new ScopedTimer("Actual clique finding");
	
	std::vector<node> pxvector(G.numberOfNodes());
	std::vector<index> pxlookup(G.upperNodeIdBound());
	uint32_t ii = 0;
	for (const node u : orderedNodes) {
		pxvector[ii] = u;
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
				pxlookup[pxvec2] = pxlookup[v];
				pxlookup[v] = xpbound - xcount - 1;
				xcount += 1;
			} else { // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
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

		std::vector<node> r = {u};
		tomita(pxvector, pxlookup, xpbound - xcount, xpbound, xpbound + pcount, r, result);

		xpbound += 1;
	}

	delete timer;

	return result;
}

void MaximalCliques::tomita(std::vector<node>& pxvector, std::vector<index>& pxlookup, uint32_t xbound, uint32_t xpbound, uint32_t pbound, std::vector<node>& r, std::vector<std::vector<node>>& result) {
	if (xbound == pbound) { //if (X, P are empty)
		result.push_back(r);
		return;
	}

	#ifndef NDEBUG
	assert(xbound >= 0);
	assert(xbound <= xpbound);
	assert(xpbound <= pbound);
	assert(pbound <= pxvector.size());
	#endif

	node u = findPivot(pxvector, pxlookup, xbound, xpbound, pbound);
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
				pxlookup[pxvec2] = pxlookup[v];
				pxlookup[v] = xpbound - xcount - 1;
				xcount += 1;
			} else if (pxlookup[v] >= xpbound && pxlookup[v] < pbound){ // v is in P
				auto pxvec2 = pxvector[xpbound + pcount];
				std::swap(pxvector[pxlookup[v]], pxvector[xpbound + pcount]);
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

		tomita(pxvector, pxlookup, xpbound - xcount, xpbound, xpbound + pcount, rplusv, result);

		auto pxvec2 = pxvector[xpbound];
		std::swap(pxvector[pxlookup[pxveci]], pxvector[xpbound]);
		pxlookup[pxvec2] = pxlookup[pxveci];
		pxlookup[pxveci] = xpbound;
		xpbound += 1;
		movedNodes.push_back(pxvector[xpbound - 1]);
	}
	
	for (node v : movedNodes) {
		//move from X -> P
		auto pxvec2 = pxvector[xpbound - 1];
		std::swap(pxvector[pxlookup[v]], pxvector[xpbound - 1]);
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
}

node MaximalCliques::findPivot(std::vector<node>& pxvector, std::vector<index>& pxlookup, uint32_t xbound, uint32_t xpbound, uint32_t pbound) {
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
