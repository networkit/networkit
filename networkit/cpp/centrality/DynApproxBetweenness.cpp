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
#include "../auxiliary/NumericTools.h"


namespace NetworKit {

DynApproxBetweenness::DynApproxBetweenness(const Graph& G, double epsilon, double delta, bool storePredecessors) : Centrality(G, true),
storePreds(storePredecessors), epsilon(epsilon), delta(delta) {

}


count DynApproxBetweenness::getNumberOfSamples() {
    return r;
}


void DynApproxBetweenness::run() {
    if (G.isDirected()) {
        throw std::runtime_error("Invalid argument: G must be undirected.");
    }
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
            sssp[i].reset(new DynDijkstra(G, u[i], storePreds));
        } else {
            sssp[i].reset(new DynBFS(G, u[i], storePreds));
        }
        DEBUG("running shortest path algorithm for node ", u[i]);
        sssp[i]->setTargetNode(v[i]);
        sssp[i]->run();
        if (sssp[i]->distances[v[i]] > 0) { // at least one path between {u, v} exists
            DEBUG("updating estimate for path ", u[i], " <-> ", v[i]);
            // random path sampling and estimation update
            sampledPaths[i].clear();
            node t = v[i];
            while (t != u[i])  {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double> > choices;
                if (storePreds) {
                    for (node z : sssp[i]->previous[t]) {
                        // workaround for integer overflow in large graphs
                        bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                        double weight;
                        tmp.ToDouble(weight);

                        choices.emplace_back(z, weight); 	// sigma_uz / sigma_us
                    }
                }
                else {
                G.forInEdgesOf(t, [&](node t, node z, edgeweight w){
                        if (Aux::NumericTools::logically_equal(sssp[i]->distances[t], sssp[i]->distances[z] + w)) {
                            // workaround for integer overflow in large graphs
                            bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                            double weight;
                            tmp.ToDouble(weight);

                            choices.emplace_back(z, weight);
                        }

                    });
                }
                INFO("Node considered: ", t);
                INFO("Source considered: ", u[i]);
                assert (choices.size() > 0);
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

    hasRun = true;

}

void DynApproxBetweenness::update(const std::vector<GraphEvent>& batch) {
    INFO ("Updating");
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
                if (storePreds) {
                    for (node z : sssp[i]->previous[t]) {
                        // workaround for integer overflow in large graphs
                        bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                        double weight;
                        tmp.ToDouble(weight);

                        choices.emplace_back(z, weight);
                    }
                }
                else {
                    G.forInEdgesOf(t, [&](node t, node z, edgeweight w){
                        if (Aux::NumericTools::logically_equal(sssp[i]->distances[t], sssp[i]->distances[z] + w)) {
                            // workaround for integer overflow in large graphs
                            bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                            double weight;
                            tmp.ToDouble(weight);

                            choices.emplace_back(z, weight);
                        }
                    });
                }
                assert (choices.size() > 0); // this should fail only if the graph is not connected
                if (choices.size() == 0) {
                    INFO ("node: ", t);
                    INFO ("source: ", u[i]);
                    INFO ("distance: ", sssp[i]->distances[t]);
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
