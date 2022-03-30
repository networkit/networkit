/*
 * GroupCloseness.hpp
 *
 *  Created on: 03.10.2016
 *      Author: elisabetta bergamini
 */

#ifndef NETWORKIT_CENTRALITY_GROUP_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_GROUP_CLOSENESS_HPP_

#include <numeric>
#include <sstream>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class GroupCloseness final : public Algorithm {
public:
    /**
     * Finds the group of nodes with highest (group) closeness centrality.
     * The algorithm is the one proposed in Bergamini et al., ALENEX 2018 and
     * finds a solution that is a (1-1/e)-approximation of the optimum.
     * The worst-case running time of this approach is quadratic, but usually
     * much faster in practice.
     *
     * @param G An unweighted graph.
     * @param k Size of the group of nodes
     * @param H If equal 0, simply runs the algorithm proposed in Bergamini et
     * al.. If > 0, interrupts all BFSs after H iterations (suggested for very
     * large networks).
     * @
     */
    GroupCloseness(const Graph &G, count k = 1, count H = 0);

    /**
     * Computes the group with maximum closeness on the graph passed in the
     * constructor.
     */
    void run() override;

    /**
     * Returns group with maximum closeness.
     */
    std::vector<node> groupMaxCloseness();

    /**
     * Computes farness (i.e., inverse of the closeness) for a given group
     * (stopping after H iterations if H > 0).
     */
    double computeFarness(const std::vector<node> &S,
                          count H = std::numeric_limits<count>::max()) const;

    /**
     * Computes the score of a specific group.
     */
    double scoreOfGroup(const std::vector<node> &group) const;

private:
    edgeweight computeImprovement(node u, count h);
    void updateDistances(node u);
    const Graph *G;
    count k = 1;
    std::vector<count> d;
    std::vector<std::vector<count>> d1Global;
    std::vector<node> S;
    count H = 0;

    void checkGroup(const std::vector<node> &group) const;
};

inline std::vector<node> GroupCloseness::groupMaxCloseness() {
    assureFinished();
    return S;
}

inline void GroupCloseness::checkGroup(const std::vector<node> &group) const {
    const count z = G->upperNodeIdBound();
    std::vector<bool> check(z, false);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(group.size()); ++i) {
        node u = group[i];
        if (u >= z) {
            std::stringstream ss;
            ss << "Error: node " << u << " is not in the graph.\n";
            throw std::runtime_error(ss.str());
        }
        if (check[u]) {
            std::stringstream ss;
            ss << "Error: the group contains duplicates of node " << u << ".\n";
            throw std::runtime_error(ss.str());
        }
        check[u] = true;
    }
}

inline double GroupCloseness::scoreOfGroup(const std::vector<node> &group) const {
    double sumDist = 0.;
    if (G->isWeighted())
        Traversal::DijkstraFrom(*G, group.begin(), group.end(),
                                [&](node, edgeweight dist) { sumDist += dist; });
    else
        Traversal::BFSfrom(*G, group.begin(), group.end(),
                           [&](node, count dist) { sumDist += static_cast<double>(dist); });

    return sumDist > 0. ? ((double)G->upperNodeIdBound() - (double)group.size()) / sumDist : 0.;
}
} /* namespace NetworKit */
#endif // NETWORKIT_CENTRALITY_GROUP_CLOSENESS_HPP_
