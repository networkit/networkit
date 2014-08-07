/*
 * DynApproxBetweenness.cpp
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

#include "DynApproxBetweenness.h"
#include "../auxiliary/Random.h"
#include "../properties/Diameter.h"
#include "../graph/Sampling.h"
#include "../graph/DynDijkstra.h"
#include "../graph/DynBFS.h"
#include "../auxiliary/Log.h"


namespace NetworKit {

DynApproxBetweenness::DynApproxBetweenness(const Graph& G, double epsilon, double delta) : Centrality(G, true), epsilon(epsilon), delta(delta) {

}

void DynApproxBetweenness::run() {
    scoreData.clear();
    scoreData.resize(G.upperNodeIdBound());
    u.clear();
    v.clear();
    sampledPaths.clear();

    double c = 0.5; // universal positive constant - see reference in paper


    edgeweight vd = Diameter::estimatedVertexDiameterPedantic(G);

    INFO("estimated diameter: ", vd);
    r = ceil((c / (epsilon * epsilon)) * (floor(log(vd - 2)) + 1 + log(1 / delta)));
    INFO("taking ", r, " path samples");
    sssp.clear();
    sssp.resize(r);
    u.resize(r);
    v.resize(r);
    sampledPaths.resize(r);

    for (count i = 0; i < r; i++) {
        DEBUG("sample ", i);
        // sample random node pair
        u[i] = Sampling::randomNode(G);
        do {
            v[i] = Sampling::randomNode(G);
        } while (v[i] == u[i]);
        if (G.isWeighted()) {
            sssp[i].reset(new DynDijkstra(G, u[i]));
        } else {
            sssp[i].reset(new DynBFS(G, u[i]));
        }
        DEBUG("running shortest path algorithm for node ", u[i]);
        sssp[i]->setTargetNode(v[i]);
        sssp[i]->run();
        if (sssp[i]->numberOfPaths(v[i]) > 0) { // at least one path between {u, v} exists
            DEBUG("updating estimate for path ", u[i], " <-> ", v[i]);
            // random path sampling and estimation update
            sampledPaths[i].clear();
            node t = v[i];
            while (t != u[i])  {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double> > choices;

                for (node z : sssp[i]->getPredecessors(t)) {
                    choices.emplace_back(z, sssp[i]->numberOfPaths(z) / (double) sssp[i]->numberOfPaths(t)); 	// sigma_uz / sigma_us
                }
                node z = Aux::Random::weightedChoice(choices);
                assert (z <= G.upperNodeIdBound());
                if (z != u[i]) {
                    scoreData[z] += 1 / (double) r;
                    sampledPaths[i].push_back(z);
                }
                t = z;
            }
        }
    }

}

void DynApproxBetweenness::update(const std::vector<GraphEvent>& batch) {
    for (node i = 0; i < r; i++) {
        sssp[i]->update(batch);
        if (sssp[i]->modified()) {
            // subtract contributions to nodes in the old sampled path
            for (node z: sampledPaths[i]) {
                scoreData[z] -= 1 / (double) r;
            }
            // sample a new shortest path
            sampledPaths[i].clear();
            node t = v[i];
            while (t != u[i])  {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double> > choices;

                for (node z : sssp[i]->getPredecessors(t)) {
                    choices.emplace_back(z, sssp[i]->numberOfPaths(z) / (double) sssp[i]->numberOfPaths(t)); 	// sigma_uz / sigma_us
                }
                node z = Aux::Random::weightedChoice(choices);
                assert (z <= G.upperNodeIdBound());
                if (z != u[i]) {
                    scoreData[z] += 1 / (double) r;
                    sampledPaths[i].push_back(z);
                }
                t = z;
            }

        }

    }
}

} /* namespace NetworKit */
