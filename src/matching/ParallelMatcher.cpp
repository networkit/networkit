/*
 * ParallelMatcher.cpp
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#include <set>

#include "ParallelMatcher.h"

namespace EnsembleClustering {

ParallelMatcher::ParallelMatcher() {
	// TODO Auto-generated constructor stub

}

ParallelMatcher::~ParallelMatcher() {
	// TODO Auto-generated destructor stub
}

Matching& ParallelMatcher::run(Graph& G) {

	// TODO: exclude isolated nodes?

	int64_t n = G.numberOfNodes();
	NodeMap<node> candidate(n, 0);					//!< candidate[v] is the preferred matching partner of v
	NodeMap<std::vector<node> > S(n, std::vector<node>()); 	//!< S[v] is a set with the potential
														//!< candidates of node v

	std::vector<node> D;	//!< targets of dominating edges
	Matching* M = new Matching(n);

	for (node v = G.firstNode(); v <= n; ++v) {

	}


	while (! D.empty()) {
		node v = D.back();
		D.pop_back();

		for (int i = 0; i < S[v].size(); ++i) {
			node x = S[v][i];
			if ((x != candidate[v]) && (! M->areMatched(x, candidate[x]))) {
				// TODO: S[x].erase()
				if (! S[x].empty()) {
					// TODO:
				}
			}
		}
	}

	return *M;
}

} /* namespace EnsembleClustering */
