/* GCE.cpp
 *
 * Created on: 06.05.2013
 * Author: cls
 */


#include "GCE.h"
#include <unordered_map>


namespace NetworKit {


GCE::GCE(const Graph& G, std::string objective) : SelectiveCommunityDetector(G), objective(objective) {
	if (G.numberOfSelfLoops() > 0) {
		throw std::runtime_error("Graphs with self-loops are not supported in GCE");
	}
}

std::map<node, std::set<node> >  GCE::run(std::set<unsigned int>& seeds) {
    std::map<node, std::set<node> > result;
    for (auto seed : seeds) {
        result[seed] = expandSeed(seed);
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
	double intWeight = 0;
    double extWeight = 0;

	struct node_property_t {
		double degInt;
		double degExt;
	};

	std::unordered_map<node, node_property_t> currentShell;

	G.forNeighborsOf(s, [&](node, node u, edgeweight ew) {
		currentShell.insert(std::make_pair(u, node_property_t {.degInt =  ew, .degExt = G.weightedDegree(u) - ew}));
		extWeight += ew;
	});


    double currentQ = 0.0; // current community quality

    auto boundary = [&](const std::set<node>& C) {
		std::set<node> sh;
		for (const auto &s : currentShell) {
			G.forNeighborsOf(s.first, [&](node u){
				if (!in(C, u)) {
					sh.insert(s.first);
				}
			});
		}
		return sh;
	};

#ifndef NDEBUG
	/**
	 * internal and external weighted degree of a node with respect to the community
	 */
	auto intExtDeg = [&](node v, const std::set<node>& C) {
		double degInt = 0;
		double degExt = 0;
		G.forNeighborsOf(v, [&](node, node u, edgeweight ew) {
			if (in(C, u)) {
				degInt += ew;
			} else {
				degExt += ew;
			}
		});
		return std::make_pair(degInt, degExt);
	};

    auto intExtWeight = [&](const std::set<node>& community) {
        double internal = 0;
        double external = 0;
        for (node u : community) {
            G.forEdgesOf(u, [&](node, node v, edgeweight ew) {
                if (in(community, v)) {
                    internal += ew;
                } else {
                    external += ew;
                }
            });
        }
        internal = internal / 2;    // internal edges were counted twice
        return std::make_pair(internal, external);
    };
#endif


    /*
     * objective function M
     * @return quality difference for the move of v to C
     */
	auto deltaM = [&](node, double degInt, double degExt, const std::set<node>& C){
		double delta = (intWeight + degInt) / (double) (extWeight - degInt + degExt);
		return delta - currentQ;
	};


    /*
     * objective function L
     * @return quality difference for the move of v to C
     */
    auto deltaL = [&](node v, double degInt, double degExt, std::set<node>& C){
    	C.insert(v);
    	double numerator = 2.0 * (intWeight + degInt) * boundary(C).size();
    	double denominator = C.size() * (extWeight - degInt + degExt);
    	C.erase(v);
        return (numerator / denominator) - currentQ;
    };

    std::function<double(node v, double degInt, double degExt, std::set<node>& C)> deltaQ;
    // select quality objective
    if (objective == "M") {
        deltaQ = deltaM;
    } else if (objective == "L") {
        deltaQ = deltaL;
    } else {
        throw std::runtime_error("unknown objective function");
    }


    // insert seed
    community.insert(s);

    // for M, quality of {s} is 0.0

	double dQMax;
	node vMax;
	do {
        // get values for current community
        assert(std::make_pair(intWeight, extWeight) == intExtWeight(community));
        // scan shell for node with maximum quality improvement
		dQMax = 0.0; 	// maximum quality improvement
		vMax = none;
//		for (node v : shell(community)) {
		for (const auto& vs : currentShell) {
            // get values for current node
			assert(intExtDeg(vs.first, community) == std::make_pair(vs.second.degInt, vs.second.degExt));

			double dQ = deltaQ(vs.first, vs.second.degInt, vs.second.degExt, community);
			TRACE("dQ: ", dQ);
			if (dQ >= dQMax) {
				vMax = vs.first;
				dQMax = dQ;
			}
		}
		TRACE("vMax: ", vMax);
		TRACE("dQMax: ", dQMax);
		if (vMax != none) {
			community.insert(vMax); 	// add best node to community

			currentShell.erase(vMax);	// remove best node from shell

			G.forNeighborsOf(vMax, [&](node, node v, edgeweight ew) { // insert external neighbors of vMax into shell
				if (!in(community, v)) {
					auto it = currentShell.find(v);
					if (it == currentShell.end()) {
						currentShell.insert(std::make_pair(v, node_property_t {.degInt = ew, .degExt = G.weightedDegree(v) - ew}));
					} else {
						it->second.degInt += ew;
						it->second.degExt -= ew;
					}

					extWeight += ew;

					assert(intExtDeg(v, community) == std::make_pair(currentShell[v].degInt, currentShell[v].degExt));
				} else {
					intWeight += ew;
					extWeight -= ew;
				}
			});
            currentQ += dQMax;     // update current community quality
			TRACE("community: ", community);
		}
	} while (vMax != none);

	return community;
}

} /* namespace NetworKit */
