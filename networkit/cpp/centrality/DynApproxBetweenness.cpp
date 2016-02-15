/*
 * DynApproxBetweenness.cpp
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

#include "DynApproxBetweenness.h"
#include "../auxiliary/Random.h"
#include "../distance/Diameter.h"
#include "../graph/Sampling.h"
#include "../graph/DynDijkstra.h"
#include "../graph/DynBFS.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/NumericTools.h"


namespace NetworKit {

DynApproxBetweenness::DynApproxBetweenness(const Graph& G, const double epsilon, const double delta, const bool storePredecessors, const double universalConstant) : Centrality(G, true),
storePreds(storePredecessors), epsilon(epsilon), delta(delta), universalConstant(universalConstant) {
  INFO("Constructing DynApproxBetweenness. storePredecessors = ", storePredecessors);
}


count DynApproxBetweenness::getNumberOfSamples() {
    return r;
}


void DynApproxBetweenness::run() {
  INFO("Inside DynApproxBetweenness. storePreds = ", storePreds);
    if (G.isDirected()) {
        throw std::runtime_error("Invalid argument: G must be undirected.");
    }
    scoreData.clear();
    scoreData.resize(G.upperNodeIdBound());
    u.clear();
    v.clear();
    sampledPaths.clear();

    Diameter diam(G, DiameterAlgo::estimatedPedantic);
    diam.run();
    edgeweight vd = diam.getDiameter().first;

    INFO("estimated diameter: ", vd);
    r = ceil((universalConstant / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 - log(delta)));
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
            INFO("Calling DynDijkstra inside run DynApproxBet");
            sssp[i].reset(new DynDijkstra(G, u[i], storePreds));
        } else {
            INFO("Calling DynBFS inside run DynApproxBet");
            sssp[i].reset(new DynBFS(G, u[i], storePreds));
        }
        DEBUG("running shortest path algorithm for node ", u[i]);

        INFO("Calling setTargetNodeon sssp instance inside run DynApproxBet");
        sssp[i]->setTargetNode(v[i]);
        INFO("Calling run on sssp instance inside run DynApproxBet");
        sssp[i]->run();
        INFO("Ran sssp");
        if (sssp[i]->distances[v[i]] > 0) { // at least one path between {u, v} exists
            DEBUG("updating estimate for path ", u[i], " <-> ", v[i]);
            INFO("Entered if statement.");
            // random path sampling and estimation update
            sampledPaths[i].clear();
            node t = v[i];
            while (t != u[i])  {
              INFO("Entered while statement");
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
                  INFO("Storepreds is false");
                  G.forInEdgesOf(t, [&](node t, node z, edgeweight w){
                        if (Aux::NumericTools::logically_equal(sssp[i]->distances[t], sssp[i]->distances[z] + w)) {
                            // workaround for integer overflow in large graphs
                            INFO("Calling number of paths");
                            bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                            INFO("Called number of paths");
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
