/*
 * Volume.cpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#include <networkit/distance/Volume.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

std::unordered_map<node, double> Volume::nodesWithinDistance(const Graph &G, double r, node n) {
    std::unordered_map<node, double> ms;
    ms.insert(std::make_pair(n, 0));
    std::vector<node> msToCheck;
    msToCheck.push_back(n);
    double r2;
    while (!msToCheck.empty()) {
        std::vector<node> msToCheckNew;
        for (auto &m : msToCheck) {
            for (node m2 : G.neighborRange(m)) {
                r2 = ms[m] + G.weight(m, m2);
                if (ms.count(m2) == 0) {
                    if (r2 <= r) {
                        ms[m2] = r2;
                        msToCheckNew.push_back(m2);
                    }
                } else {
                    ms[m2] = std::fmin(ms[m2], r2);
                }
            }
        }
        msToCheck = msToCheckNew;
    }
    return ms;
}

double Volume::volume(const Graph &G, double r, count samples) {
    double x = 0;
    for (count j = 0; j < samples; j++) {
        x += Volume::nodesWithinDistance(G, r, GraphTools::randomNode(G)).size();
    }
    return x / samples;
}

std::vector<double> Volume::volume(const Graph &G, std::vector<double> rs, count samples) {
    std::vector<double> xs(rs.size(), 0);
    double rmax = *std::max_element(std::begin(rs), std::end(rs));
    for (count j = 0; j < samples; j++) {
        std::unordered_map<node, double> ms =
            Volume::nodesWithinDistance(G, rmax, GraphTools::randomNode(G));
        for (count i = 0; i < rs.size(); ++i)
            xs[i] += std::count_if(ms.begin(), ms.end(), [&](const auto &nodeDist) -> bool {
                return nodeDist.second <= rs[i];
            });
    }
    std::vector<double> ys;
    ys.reserve(xs.size());
    std::transform(xs.cbegin(), xs.cend(), std::back_inserter(ys),
                   [samples](double x) { return x / static_cast<double>(samples); });
    return ys;
}

} /* namespace NetworKit */
