/* GCE.cpp
 *
 * Created on: 06.05.2013
 * Author: cls
 */


#include "GCE.h"


namespace NetworKit {


GCE::GCE(const Graph& G) : SelectiveCommunityDetector(G) {

}

void GCE::run(std::set<unsigned int>& seeds) {

}

std::set<node> GCE::expandSeed(node s) {
	auto in = [](std::set<node> A, node x) {
		return (A.find(x) != A.end());
	};

	std::set<node> community;
	community.insert(s);

	// values per community
	count intEdges = 0;
	count extEdges = 0;

	/** @return the shell of the given community */
	auto shell = [&](std::set<node> C) {
		std::set<node> sh;
		for (node v : C) {
			G.forNeighborsOf(v, [&](node u){
				if (!in(C, u)) {
					sh.insert(u);
				}
			});
		}
		return sh;
	};

	/** 
	 * internal and external degree of a node with respect to the community
	 */
	auto intExtDeg = [&](node v, std::set<node> C) {
		count degInt = 0;
		count degExt = 0;
		G.forNeighborsOf(v, [&](node u) {
			if (in(C, u)) {
				degInt += 1;
			} else {
				degExt += 1;
			}
		});
		return std::make_pair(degInt, degExt);
	};


	auto deltaM = [&](node v, std::set<node> C){
		count degInt, degExt;
		std::tie(degInt, degExt) =  intExtDeg(v, C);
		double delta = (intEdges + degInt) / (double) (extEdges + degInt + degExt);
		return delta;
	};


	auto deltaQ = deltaM; // select quality objective


	double dQMax;
	node vMax;
	do {
		// scan shell for node with maximum quality improvement
		dQMax = 0.0; 	// maximum quality improvement
		vMax = none;
		for (node v : shell(community)) {
			double dQ = deltaQ(v, community);
			TRACE("dQ: ", dQ);
			if (dQ >= dQMax) {
				vMax = v;
				dQMax = dQ;
			}
		}
		TRACE("vMax: ", vMax);
		TRACE("dQMax: ", dQMax);
		if (vMax != none) {
			community.insert(vMax); 	// add best node to community
			TRACE("community: ", community);
		}
	} while (vMax != none);

	return community;
}

} /* namespace NetworKit */

