/*
 * DegreeSequenceGenerator_Hoske.cpp
 *
 *  Created on: 05-12-2013
 *      Author: dhoske
 */

#include "DegreeSequenceGenerator_Hoske.h"

namespace NetworKit {

Graph DegreeSequenceGenerator_Hoske::generate() {
	using namespace std;
    Graph G(n);

    /* Pairs of (remaining degree, node id). */
    typedef pair<count, node> DegNode;
    vector<DegNode> rem_deg(n);
    for (node v = 0; v < n; ++v) {
    	rem_deg[v] = {deg_seq[v], v};
    }

    /* In each loop: fulfil degree requirement for one node. */
    for (count k = 0; k < n; ++k) {
    	/* Distribute edges for node with largest remaining degree. */
    	sort(begin(rem_deg) + k, end(rem_deg), greater<DegNode>());
    	count& deg = rem_deg[k].first;
    	node v = rem_deg[k].second;
    	if (deg > n - k - 1) {
    		throw invalid_argument("Degree sequence not realizable.");
    	}

    	deg = 0;
    	for (count l = k + 1; l <= k + deg; ++l) {
    		rem_deg[l].first--;
    		G.addEdge(v, rem_deg[l].second);
    	}
    }

    return G;
}

} /* namespace NetworKit */

