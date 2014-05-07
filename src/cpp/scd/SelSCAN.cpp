#include "SelSCAN.h"

#include <queue>

namespace NetworKit {

void SelSCAN::run(std::set<unsigned int>& seeds) {

	auto dist = [&](node u, node v) {
		// TODO: implement
		return 0.0;
	};

	/** 
	 * @return the epsilon-neighborhood of a node
	 */
	auto epsilonNeighborhood = [&](node u) {
		std::set<node> N;
		G.forNeighborsOf(u, [&](node v) {
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
		return (epsilonNeighborhood(u).size() >= kappa);
	};


	auto coreSearch = [&](std::queue<node>& Q) {
		while (!Q.empty()) {
			node x = Q.front(); Q.pop();
			// TODO: assign x to current community
			if (core(x)) {
				for (node y : epsilonNeighborhood(x)) {
					// TODO: loop body
				}
			}
		}

	};

	for (auto s: seeds) {
		std::queue<node> Q;
		if (core(s)) {
			// TODO: assign s to new community
			Q.push(s);
			coreSearch(Q);
		} else {
			// find core with minimum distance to s in epsilon neighborhood of s
			node minC = none;
			double minD = std::numeric_limits<double>::max();
			for (node c : epsilonNeighborhood(s)) {
				if (core(c)) {
					if (dist(s, c) < minD) {
						minD = dist(s, c);
						minC = c;
					}
				}
			}

			if (minC == none) { 	// no core in neighborhood
				// TODO: assign s as outlier
			} else {
				// TODO: assign s to new community
				Q.push(s);
				coreSearch(Q);
			}
		} // 
	}

}


}