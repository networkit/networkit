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
	NodeMap<std::set<node> > S(n, std::set<node>()); 	//!< S[v] is a set with the potential
																		//!< candidates of node v

	std::set<node> D;	//!< targets of dominating edges
	Matching M(n);

	// TODO: for all nodes

	return M;
}

} /* namespace EnsembleClustering */
