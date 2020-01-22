/*
 * PageRankNibble.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include <algorithm>
#include <unordered_set>
#include <vector>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/scd/ApproximatePageRank.hpp>
#include <networkit/scd/PageRankNibble.hpp>

namespace NetworKit {

PageRankNibble::PageRankNibble(const Graph& g, double alpha, double epsilon): SelectiveCommunityDetector(g), alpha(alpha), epsilon(epsilon) {}

std::set<node> PageRankNibble::bestSweepSet(std::vector<std::pair<node, double>>& pr) {
    TRACE("Finding best sweep set. Support size: ",  pr.size());

    // order vertices
    TRACE("Before sorting");
    for (size_t i = 0; i < pr.size(); i++) {
        pr[i].second = pr[i].second / G->weightedDegree(pr[i].first, true);
    }
    auto comp([&](const std::pair<node, double>& a, const std::pair<node, double>& b) {
        return a.second > b.second;
    });
    Aux::Parallel::sort(pr.begin(), pr.end(), comp);
    TRACE("After sorting");

    #ifndef NDEBUG
    for (auto it = pr.begin(); it != pr.end(); it++) {
        TRACE("(", it->first, ", ", it->second, ")");
    }
    #endif

    // find best sweep set w.r.t. conductance
    double bestCond = std::numeric_limits<double>::max();
    double cut = 0.0;
    double volume = 0.0;
    index bestSweepSetIndex = 0;
    std::unordered_set<node> withinSweepSet;
    std::vector<node> currentSweepSet;

    // generate total volume.
    double totalVolume = G->totalEdgeWeight() * 2;

    for (auto it = pr.begin(); it != pr.end(); it++) {
        // update sweep set
        node v = it->first;
        double wDegree = 0.0;
        G->forNeighborsOf(v, [&](node, node neigh, edgeweight w) {
            wDegree += w;
            if (withinSweepSet.find(neigh) == withinSweepSet.end()) {
                cut += w;
            } else {
                cut -= w;
            }
        });
        volume += wDegree;
        currentSweepSet.push_back(v);
        withinSweepSet.insert(v);

        // compute conductance
        double cond = cut / std::min(volume, totalVolume - volume);

        if ((cond < bestCond) && (currentSweepSet.size() < G->numberOfNodes())) {
            bestCond = cond;
            bestSweepSetIndex = currentSweepSet.size();
        }
    }

    DEBUG("Best conductance: ", bestCond, "\n");

    std::set<node> bestSweepSet(currentSweepSet.begin(), currentSweepSet.begin() + bestSweepSetIndex);
    return bestSweepSet;
}

std::set<node> PageRankNibble::expandSeed(node seed) {
    DEBUG("APR(G, ", alpha, ", ", epsilon, ")");
    ApproximatePageRank apr(*G, alpha, epsilon);
    std::vector<std::pair<node, double>> pr = apr.run(seed);
    return bestSweepSet(pr);
}

std::map<node, std::set<node> >  PageRankNibble::run(const std::set<node>& seeds) {
    std::map<node, std::set<node> > result;
    for (auto seed : seeds) {
        auto community = expandSeed(seed);
        result[seed] = community;
    }
    return result;
}

} /* namespace NetworKit */
