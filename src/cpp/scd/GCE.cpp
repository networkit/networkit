/* GCE.cpp
 *
 * Created on: 06.05.2013
 * Author: cls
 */


#include "GCE.h"


namespace NetworKit {


GCE::GCE(const Graph& G) : SelectiveCommunityDetector(G) {

}

std::map<node, std::set<node> >  GCE::run(std::set<unsigned int>& seeds) {
    std::map<node, std::set<node> > result;
    for (auto seed : seeds) {
        auto community = expandSeed(seed);
        result[seed] = community;
    }
    return result;
}

std::set<node> GCE::expandSeed(node s) {
    /**
    * Check if set contains node.
    */
	auto in = [](const std::set<node>& A, node x) {
		return (A.find(x) != A.end());
	};

	std::set<node> community;

	// values per community
	count intEdges = 0;
    count extEdges = 0;

    count bSize = 0; // boundary size
    count dBSize = 0; // TODO:

    double currentQ = 0.0; // current community quality

    // values per node
    count degInt, degExt; // internal, external degree



	/** @return the shell of the given community */
	auto shell = [&](const std::set<node>& C) {
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
	auto intExtDeg = [&](node v, const std::set<node>& C) {
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


    auto intExtEdges = [&](const std::set<node>& community) {
        count internal = 0;
        count external = 0;
        for (node u : community) {
            G.forEdgesOf(u, [&](node u, node v) {
                if (in(community, v)) {
                    internal += 1;
                } else {
                    external += 1;
                }
            });
        }
        internal = internal / 2;    // internal edges were counted twice
        return std::make_pair(internal, external);
    };


    /*
     * objective function M
     * @return quality difference for the move of v to C
     */
	auto deltaM = [&](node v, const std::set<node>& C){
		double delta = (intEdges + degInt) / (double) (extEdges - degInt + degExt);
		return delta - currentQ;
	};


    /*
     * objective function L
     * @return quality difference for the move of v to C
     */
    auto deltaL = [&](node v, const std::set<node>& C){
        double numerator = (extEdges + 2 * degInt) / (double) (community.size() + 1);
        double denominator = (extEdges - degInt + degExt) / (bSize + dBSize);
        return (numerator / denominator) - currentQ;
    };



	auto deltaQ = deltaM; // select quality objective


    // insert seed
    community.insert(s);

    // for M, quality of {s} is 0.0

	double dQMax;
	node vMax;
	do {
        // get values for current community
        std::tie(intEdges, extEdges) = intExtEdges(community);
        // scan shell for node with maximum quality improvement
		dQMax = 0.0; 	// maximum quality improvement
		vMax = none;
		for (node v : shell(community)) {
            // get values for current node
            std::tie(degInt, degExt) =  intExtDeg(v, community);
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
            currentQ += dQMax;     // update current community quality
			TRACE("community: ", community);
		}
	} while (vMax != none);

	return community;
}

} /* namespace NetworKit */
