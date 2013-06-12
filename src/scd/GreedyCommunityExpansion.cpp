/*
 * GreedyCommunityExpansion.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GreedyCommunityExpansion.h"

namespace NetworKit {

GreedyCommunityExpansion::GreedyCommunityExpansion() {
	// TODO Auto-generated constructor stub

}

GreedyCommunityExpansion::~GreedyCommunityExpansion() {
	// TODO Auto-generated destructor stub
}

std::unordered_set<node> GreedyCommunityExpansion::run(Graph& G, node s) {


	std::unordered_set<node> community;
	community.insert(s); // begin with C_s = {s}

	std::unordered_set<node> shell; // shell are the nodes outside of the community with edges to nodes inside
	// initialize shell to N(s)
	G.forNeighborsOf(s, [&](node v) {
		shell.insert(v);
	});

	double deltaQMax;
	node vMax;

	do {
		for (node v : shell) { // TODO: optionally order the nodes by acceptability
			// TODO: evaluate Delta Q

		}
		// TODO: if Delta Q* > 0
		community.insert(vMax);
		// update shell incrementally in O(deg(vMax))
		shell.erase(vMax);
		G.forNeighborsOf(vMax, [&](node v){
			if (community.find(v) == community.end()) {
				// v is not in C
				shell.insert(v);
			}
		});
	} while (deltaQMax > 0);


	// TODO: optional trimming phase according to node fitness
	return community;
}

} /* namespace NetworKit */
