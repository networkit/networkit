#include "SelSCAN.h"

#include <queue>

namespace NetworKit {

SelSCAN::SelSCAN(const Graph& G, count kappa, double epsilon) : SelectiveCommunityDetector(G), kappa(kappa), epsilon(epsilon), algebraicDistance(NULL), intersector(G.upperNodeIdBound()) {
	/**
	 * expresses neighborhood overlap
	 */
	auto neighborhoodOverlap = [&](node u, node v) -> double {
		auto Nu = G.neighbors(u);
		auto Nv = G.neighbors(v);
		Nu.push_back(u);
		Nv.push_back(v);
		// std::set<node> intersect;
		// set_intersection(Nu.begin(), Nu.end(), Nv.begin(), Nv.end(), std::inserter(intersect, intersect.begin()));
		auto intersect = intersector.intersect(Nu, Nv);

		return 1 - (intersect.size() / sqrt(Nu.size() * Nv.size()));
	};

	dist = neighborhoodOverlap;

}

SelSCAN::SelSCAN(const Graph& G, count kappa, double epsilon, AlgebraicDistance& ad) : SelectiveCommunityDetector(G), kappa(kappa), epsilon(epsilon), algebraicDistance(&ad), intersector(G.upperNodeIdBound()) {
	auto algDist = [&](node u, node v) -> double {
		return algebraicDistance->distance(u, v);
	};

	dist = algDist;
}

std::map<node, std::set<node> >  SelSCAN::run(std::set<unsigned int>& seeds) {

	/**
	 * @return the epsilon-neighborhood of a node
	 */
	auto epsilonNeighborhood = [&](node u) {
		std::set<node> N;
		G.forNeighborsOf(u, [&](node v) {
			TRACE("neighbor ", v, " has distance ", dist(u, v));
			if (dist(u, v) < epsilon) {
				N.insert(v);
			}
		});
		return N;
	};

	/**
	 * determines if a node is a core
	 */
	auto core = [&](node u) {
		TRACE(u, " has an epsilon neighborhood of size ", epsilonNeighborhood(u).size() );
		return (epsilonNeighborhood(u).size() >= kappa);
	};





	// map from node to community label
	std::unordered_map<node, index> eta;
	index label = 0;	// current community label

	auto newLabel = [&]() {
		label += 1;
		return label;
	};


	/**
	* breadth-first search which appends neighbors to the community
	* and cores to the queue.
	*/
	auto coreSearch = [&](std::queue<node>& Q, index label) {
		TRACE("starting core search");
		while (!Q.empty()) {
			node x = Q.front(); Q.pop();
			assert (core(x));
			for (node y : epsilonNeighborhood(x)) {
				// assert ((eta.find(y) == eta.end()) || (eta[y] == none));	// not assigned or outlier - FAILS
				if (eta.find(y) == eta.end()) {
					if (core(y)) {
						Q.push(y);
					}
					eta[y] = label; // add to community
					TRACE(y, " added to community ", eta[y]);
				}
				// labelled neighbors are ignored
			}
		}

	};


	for (auto s: seeds) {
		DEBUG("seed ", s);
		if (eta.find(s) == eta.end()) {	// s not yet assigned to community
			std::queue<node> Q;
			if (core(s)) {	// if s is a core, start new community and grow ith by core search
				TRACE(s, " is a core");
				eta[s] = newLabel();
				TRACE(s, " added to community ", eta[s]);
				Q.push(s);
				coreSearch(Q, eta[s]);
			} else {  // if s is not a core, find core with minimum distance to s in epsilon neighborhood of s
				TRACE(s, " is not a core");
				node minC = none;
				double minD = std::numeric_limits<double>::max();
				for (node c : epsilonNeighborhood(s)) {	// look for closest core
					if (core(c)) {
						if (dist(s, c) < minD) {
							minD = dist(s, c);
							minC = c;
						}
					}
				}

				if (minC == none) { 	// no core in neighborhood
					TRACE(s, " is an outlier");
					eta[s] = none; 	// mark s as outlier
				} else {		// start new community for s and grow it by core search from c
					TRACE(s, " has closest core neighbor ", minC);
					eta[s] = newLabel();
					TRACE(s, " added to community ", eta[s]);
					Q.push(minC);
					coreSearch(Q, eta[s]);
				}
			} //
		}

	} // end community discovery

	DEBUG("eta: ", eta);

	// extract sets from labels
	std::map<node, std::set<node> > result;
	for (auto s : seeds) {
		index ls = eta[s];
		std::set<node> community;
		if (ls != none) {
			for (auto kv : eta) {
				if (kv.second == ls) {
					community.insert(kv.first);
				}
			}
		}
		result[s] = community;
	}
	return result;

}


}
