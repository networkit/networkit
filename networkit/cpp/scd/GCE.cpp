/* GCE.cpp
 *
 * Created on: 06.05.2013
 * Author: cls
 */

#include <unordered_map>

#include <networkit/scd/GCE.hpp>

namespace {
    template <bool> struct ConditionalCount {
        void increaseNumBoundaryNeighbors() {}
        void decreaseNumBoundaryNeighbors() {}
        NetworKit::count getNumBoundaryNeighbors() {
            throw std::logic_error("Error, num boundary neighbors not available in false branch");
        }
    };
    template <> struct ConditionalCount<true> {
        NetworKit::count numBoundaryNeighbors;
        void increaseNumBoundaryNeighbors() {
            numBoundaryNeighbors += 1;
        }
        void decreaseNumBoundaryNeighbors() {
            assert(numBoundaryNeighbors > 0);
            numBoundaryNeighbors -= 1;
        }
        NetworKit::count getNumBoundaryNeighbors() {
            return numBoundaryNeighbors;
        }
    };
}

namespace NetworKit {


GCE::GCE(const Graph& G, std::string objective) : SelectiveCommunityDetector(G), objective(objective) {
    if (G.numberOfSelfLoops()) {
        throw std::runtime_error("Graphs with self-loops are not supported in GCE");
    }
}

std::map<node, std::set<node> >  GCE::run(const std::set<node>& seeds) {
    std::map<node, std::set<node> > result;
    for (auto seed : seeds) {
        result[seed] = expandSeed(seed);
    }
    return result;
}


template <bool objectiveIsM>
std::set<node> expandseed_internal(const Graph&G, node s) {
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

    struct node_property_t : ConditionalCount<!objectiveIsM> {
        double degInt;
        double degExt;
    };

    std::unordered_map<node, node_property_t> currentShell;

    double currentQ = 0.0; // current community quality

    // Current boundary. Stores for every node in the boundary how many neighbors it has outside the community.
    std::unordered_map<node, count> currentBoundary;

#ifndef NDEBUG
    // The boundary is defined as all nodes of C that have a neighbor not in C
    auto boundary = [&](const std::set<node>& C) {
        std::set<node> sh;
        for (node u : C) {
            G.forNeighborsOf(u, [&](node v){
                if (!in(C, v)) {
                    sh.insert(u);
                }
            });
        }
        return sh;
    };

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
        internal = internal / 2;	// internal edges were counted twice
        return std::make_pair(internal, external);
    };
#endif

    auto addNodeToCommunity = [&](node u) {
        community.insert(u); 	// add node to community

        currentShell.erase(u);	// remove node from shell

        node boundaryNeighbor = none; // for L: if u is in the boundary and has only one neighbor outside of the community, store it here.
        auto boundaryIt = currentBoundary.end();

        G.forNeighborsOf(u, [&](node, node v, edgeweight ew) { // insert external neighbors of u into shell
            if (!in(community, v)) {
                auto it = currentShell.find(v);
                if (it == currentShell.end()) {
                    std::tie(it, std::ignore) = currentShell.insert({v, node_property_t {}});
                    it->second.degExt = G.weightedDegree(v);
                }

                it->second.degInt += ew;
                it->second.degExt -= ew;

                extWeight += ew;
                if (!objectiveIsM) {
                    if (boundaryIt == currentBoundary.end()) {
                        std::tie(boundaryIt, std::ignore) = currentBoundary.insert({u, 0});
                        boundaryNeighbor = v;
                    }

                    ++boundaryIt->second;
                }

                assert(intExtDeg(v, community) == std::make_pair(currentShell[v].degInt, currentShell[v].degExt));
            } else {
                if (!objectiveIsM) {
                    auto it = currentBoundary.find(v);
                    assert(it != currentBoundary.end());
                    it->second -= 1;
                    if (it->second == 0) {
                        currentBoundary.erase(it);
                    } else if (it->second == 1) {
                        G.forNeighborsOf(v, [&](node x) {
                            auto it = currentShell.find(x);
                            if (it != currentShell.end()) {
                                it->second.increaseNumBoundaryNeighbors();
                            }
                        });
                    }
                }

                intWeight += ew;
                extWeight -= ew;
            }
        });

        if (!objectiveIsM && boundaryIt != currentBoundary.end() && boundaryIt->second == 1) {
            assert(boundaryNeighbor != none);
            currentShell[boundaryNeighbor].increaseNumBoundaryNeighbors();
        }

        assert(objectiveIsM || boundary(community).size() == currentBoundary.size());
    };

    addNodeToCommunity(s);

    /*
     * objective function M
     * @return quality difference for the move of v to C
     */
    auto deltaM = [&](node, double degInt, double degExt, const std::set<node> &){
        double delta = (intWeight + degInt) / (double) (extWeight - degInt + degExt);
        return delta - currentQ;
    };


    /*
     * objective function L
     * @return quality difference for the move of v to C
     */
    auto deltaL = [&](node v, double degInt, double degExt, std::set<node>& C){
    // Compute difference in boundary size: for each neighbor where we are the last
    // external neighbor decrease by 1, if v has an external neighbor increase by 1
    int64_t boundary_diff = 0;
    if (degExt > 0) {
        boundary_diff += 1;
    }

    boundary_diff -= currentShell[v].getNumBoundaryNeighbors();

#ifndef NDEBUG
    int64_t boundary_diff_debug = 0;
    bool v_in_boundary = false;
    G.forNeighborsOf(v, [&](node x) {
        auto it = currentBoundary.find(x);
        if (it != currentBoundary.end()) {
            if (it->second == 1) {
                boundary_diff_debug -= 1;
            }
        } else if (!v_in_boundary) {
            boundary_diff_debug += 1;
            v_in_boundary = true;
        }
    });

    assert(boundary_diff == boundary_diff_debug);
#endif
    double numerator = 2.0 * (intWeight + degInt) * (currentBoundary.size() + boundary_diff);
    double denominator = (C.size() + 1) * (extWeight - degInt + degExt);
        return (numerator / denominator) - currentQ;
    };

    // select quality objective
    auto deltaQ = [&](node v, double degInt, double degExt, std::set<node>& C) -> double {
        if (objectiveIsM) {
            return deltaM(v, degInt, degExt, C);
        } else {
            return deltaL(v, degInt, degExt, C);
        }
    };

    // for M, quality of {s} is 0.0

    double dQMax;
    node vMax;
    do {
        // get values for current community
        assert(std::make_pair(intWeight, extWeight) == intExtWeight(community));
        // scan shell for node with maximum quality improvement
        dQMax = 0.0; 	// maximum quality improvement
        vMax = none;
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
            addNodeToCommunity(vMax);  // add best node to community
            currentQ += dQMax;   // update current community quality
            TRACE("community: ", community);
        }
    } while (vMax != none);

    return community;
}

std::set<node> GCE::expandSeed(node s) {
    if (objective == "M") {
        return expandseed_internal<true>(*G, s);
    } else if (objective == "L") {
        return expandseed_internal<false>(*G, s);
    } else {
        throw std::runtime_error("unknown objective function");
    }
}

} /* namespace NetworKit */
