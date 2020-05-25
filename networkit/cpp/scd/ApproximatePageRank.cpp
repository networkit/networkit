/*
 * ApproximatePageRank.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include <queue>

#include <networkit/scd/ApproximatePageRank.hpp>

namespace NetworKit {

ApproximatePageRank::ApproximatePageRank(const Graph& g, double alpha_, double epsilon):
        G(&g), alpha(alpha_), eps(epsilon) {}

std::vector<std::pair<node, double>> ApproximatePageRank::run(node seed) {
    pr_res[seed] = std::make_pair(0.0, 1.0);
    std::queue<node> activeNodes;
    activeNodes.push(seed);

    auto push = [&](const node u, std::queue<node> &activeNodes) {
        double res = pr_res[u].second;
        double volume = G->weightedDegree(u, true);

        G->forNeighborsOf(u, [&](node, const node v, const edgeweight w) {
            double mass = (1.0 - alpha) * res * w / (2.0 * volume);
            double vol_v = G->weightedDegree(v, true);
            // the first check is for making sure the node is not added twice.
            // the second check ensures that enough residual is left.
            if ( pr_res[v].second < vol_v * eps && (pr_res[v].second + mass) >= eps * vol_v ) {
                activeNodes.push(v);
            }
            pr_res[v].second += mass;
        });

        pr_res[u] = std::make_pair(pr_res[u].first + alpha * res, (1.0 - alpha) * res / 2);
        if ((pr_res[u].second / volume) >= eps) {
            activeNodes.push(u);
        }
    };

    while (!activeNodes.empty()) {
        node v =  activeNodes.front();
        activeNodes.pop();
        push(v, activeNodes);
    }

    std::vector<std::pair<node, double>> pr;
    pr.reserve(pr_res.size());

    for (auto it = pr_res.begin(); it != pr_res.end(); it++) {
        pr.push_back(std::make_pair(it->first, it->second.first));
    }

    return pr;
}

} /* namespace NetworKit */
