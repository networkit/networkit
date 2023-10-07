/*
 * DynApproxBetweenness.cpp
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/DynBFS.hpp>
#include <networkit/distance/DynDijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/TopologicalSort.hpp>

namespace NetworKit {

DynApproxBetweenness::DynApproxBetweenness(const Graph &G, double epsilon, double delta,
                                           bool storePredecessors, double universalConstant)
    : Centrality(G, true), epsilon(epsilon), delta(delta), storePreds(storePredecessors),
      universalConstant(universalConstant) {}

std::vector<node>
DynApproxBetweenness::sortComponentsTopologically(Graph &sccDAG, StronglyConnectedComponents &scc) {
    G.forEdges([&](node u, node v) {
        if (scc.componentOfNode(u) != scc.componentOfNode(v)) {
            if (!sccDAG.hasEdge(scc.componentOfNode(u), scc.componentOfNode(v))) {
                sccDAG.addEdge(scc.componentOfNode(u), scc.componentOfNode(v));
            }
        }
    });

    TopologicalSort topSort(sccDAG);
    topSort.run();
    return topSort.getResult();
}

count DynApproxBetweenness::computeVDdirected() {
    StronglyConnectedComponents scc(G);
    scc.run();
    count sccNum = scc.numberOfComponents();
    Graph sccDAG(sccNum, false, true);

    std::vector<node> sorted = sortComponentsTopologically(sccDAG, scc);

    std::vector<bool> marked1(G.upperNodeIdBound()), marked2(G.upperNodeIdBound());
    std::vector<node> sources(sccNum, 0);
    G.forNodes([&](node u) {
        count x = scc.componentOfNode(u);
        if (sources[x] == 0) {
            sources[x] = u;
        }
    });
    // we compute the bound on VD for each SCC
    std::vector<count> vdapproxs(sccNum);
    std::queue<node> q, qNext;
    for (count i = 0; i < sccNum; i++) {
        count source = sources[sorted[i]];
        // forward BFS
        marked1[source] = true;
        count dist1 = 0;
        q.push(source);
        do {
            node u = q.front();
            q.pop();
            G.forNeighborsOf(u, [&](node v) {
                if (!marked1[v] && scc.componentOfNode(u) == scc.componentOfNode(v)) {
                    qNext.push(v);
                    marked1[v] = true;
                }
            });
            if (q.empty() && !qNext.empty()) {
                q.swap(qNext);
                ++dist1;
            }
        } while (!q.empty());

        // backward BFS
        marked2[source] = true;
        count dist2 = 0;
        q.push(source);
        do {
            node u = q.front();
            q.pop();
            G.forInNeighborsOf(u, [&](node v) {
                if (!marked2[v] && scc.componentOfNode(u) == scc.componentOfNode(v)) {
                    qNext.push(v);
                    marked2[v] = true;
                }
            });
            if (q.empty() && !qNext.empty()) {
                q.swap(qNext);
                ++dist2;
            }
        } while (!q.empty());
        // source node is not taken into account in any of the two approximations
        vdapproxs[sorted[i]] = dist1 + dist2 + 1;
        if (vdapproxs[sorted[i]] == 0) // for singletons
            vdapproxs[sorted[i]] = 1;
    }
    // starting from the bottom, we compute all the cumulative vd approximations
    count vd = vdapproxs[0];
    std::vector<count> vdSuccessors(sccNum, 0);
    for (count i = sccNum; i > 0; i--) {
        node c = sorted[i - 1];
        if (sccDAG.degreeOut(c) > 0) {
            vdapproxs[c] += vdSuccessors[c];
        }
        if (vdapproxs[c] > vd) {
            vd = vdapproxs[c];
        }
        sccDAG.forInEdgesOf(c, [&](node c_pred) {
            if (vdSuccessors[c_pred] < vdapproxs[c])
                vdSuccessors[c_pred] = vdapproxs[c];
        });
    }
    return vd;
}

count DynApproxBetweenness::getNumberOfSamples() const noexcept {
    return r;
}

void DynApproxBetweenness::run() {
    scoreData.clear();
    scoreData.resize(G.upperNodeIdBound());
    u.clear();
    v.clear();
    sampledPaths.clear();
    count vd;
    if (G.isDirected()) {
        vd = computeVDdirected();
    } else {
        Diameter diam(G, DiameterAlgo::ESTIMATED_PEDANTIC);
        diam.run();
        vd = diam.getDiameter().first;
    }
    r = std::ceil((universalConstant / (epsilon * epsilon))
                  * (tlx::integer_log2_floor(vd - 2) + 1.0 + std::log(1.0 / delta)));
    sssp.clear();
    sssp.resize(r);
    sampleNewPaths(0, r);
    hasRun = true;
}

void DynApproxBetweenness::update(GraphEvent e) {
    std::vector<GraphEvent> batch = {e};
    updateBatch(batch);
}

void DynApproxBetweenness::updateBatch(const std::vector<GraphEvent> &batch) {
    assert(sssp.size() == r);
    for (node i = 0; i < r; i++) {
        sssp[i]->updateBatch(batch);
        if (!sssp[i]->modified()) // Skip unaffected node i
            continue;
        // subtract contributions to nodes in the old sampled path
        for (node z : sampledPaths[i]) {
            scoreData[z] -= 1 / (double)r;
        }
        if (sssp[i]->distance(v[i]) == infDist) // Skip unreachable node v[i]
            continue;
        // sample a new shortest path
        sampledPaths[i].clear();
        node t = v[i];
        while (t != u[i]) {
            // sample z in P_u(t) with probability sigma_uz / sigma_us
            std::vector<std::pair<node, double>> choices;
            G.forInEdgesOf(t, [&](node t, node z, edgeweight w) {
                if (Aux::NumericTools::logically_equal(sssp[i]->distance(t),
                                                       sssp[i]->distance(z) + w)) {
                    // workaround for integer overflow in large graphs
                    double weight =
                        (sssp[i]->getNumberOfPaths(z) / sssp[i]->getNumberOfPaths(t)).ToDouble();
                    choices.emplace_back(z, weight);
                }
            });
            assert(!choices.empty()); // this should fail only if the graph is not connected
            node z = Aux::Random::weightedChoice(choices);
            assert(G.hasNode(z));
            if (z != u[i]) {
                scoreData[z] += 1 / (double)r;
                sampledPaths[i].push_back(z);
            }
            t = z;
        }
    }
    count new_vd;
    if (G.isDirected()) {
        new_vd = computeVDdirected();
    } else { // case for undirected already present
        Diameter diam(G, DiameterAlgo::ESTIMATED_PEDANTIC);
        diam.run();
        new_vd = diam.getDiameter().first;
    }
    count new_r = std::ceil((universalConstant / (epsilon * epsilon))
                            * (tlx::integer_log2_floor(new_vd - 2) + 1.0 + std::log(1.0 / delta)));
    if (new_r > r) {
        sampleNewPaths(r, new_r);
        // multiply all scores by r/new_r
        G.forNodes([&](node n) { scoreData[n] *= r / double(new_r); });
        r = new_r;
    }
}

void DynApproxBetweenness::sampleNewPaths(count start, count end) {
    for (count i = start; i < end; i++) {
        // sample random node pair
        node u1, u2;
        u1 = GraphTools::randomNode(G);
        do {
            u2 = GraphTools::randomNode(G);
        } while (u1 == u2);
        u.push_back(u1);
        v.push_back(u2);
        if (G.isWeighted()) {
            sssp[i] = std::make_unique<DynDijkstra>(G, u[i], storePreds);
        } else {
            sssp[i] = std::make_unique<DynBFS>(G, u[i], storePreds);
        }
        sssp[i]->run();
        std::vector<node> path;
        if (sssp[i]->distance(v[i]) < infDist) { // at least one path between {u, v} exists
            // random path sampling and estimation update
            node t = v[i];
            while (t != u[i]) {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double>> choices;
                if (storePreds) {
                    for (node z : sssp[i]->previous[t]) {
                        // workaround for integer overflow in large graphs
                        double weight =
                            (sssp[i]->getNumberOfPaths(z) / sssp[i]->getNumberOfPaths(t))
                                .ToDouble();
                        assert(weight != 0);
                        choices.emplace_back(z, weight); // sigma_uz / sigma_us
                    }
                } else {
                    G.forInEdgesOf(t, [&](node t, node z, edgeweight w) {
                        if (Aux::NumericTools::logically_equal(sssp[i]->distance(t),
                                                               sssp[i]->distance(z) + w)) {
                            // workaround for integer overflow in large graphs
                            double weight =
                                (sssp[i]->getNumberOfPaths(z) / sssp[i]->getNumberOfPaths(t))
                                    .ToDouble();
                            assert(weight != 0);
                            choices.emplace_back(z, weight);
                        }
                    });
                }

                assert(!choices.empty());
                node z = Aux::Random::weightedChoice(choices);
                assert(G.hasNode(z));
                if (z != u[i]) {
                    scoreData[z] += 1 / (double)r;
                    path.push_back(z);
                }
                t = z;
            }
        }
        sampledPaths.emplace_back(std::move(path));
    }
}

constexpr edgeweight DynApproxBetweenness::infDist;

} // namespace NetworKit
