/*
 * DynBetweenness.cpp
 *
 *  Created on: 29.07.2014
 *      Author: ebergamini
 */
#include <stack>
#include <queue>
#include <memory>

#include "../auxiliary/PrioQueue.h"
#include "DynBetweenness.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/NumericTools.h"
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {
DynBetweenness::DynBetweenness(const Graph& G, bool storePredecessors) : Centrality(G, normalized),
maxDistance(G.upperNodeIdBound()),
npaths(G.upperNodeIdBound(), std::vector<bigfloat>(G.upperNodeIdBound())),
distances(G.upperNodeIdBound(), std::vector<edgeweight>(G.upperNodeIdBound())),
dependencies(G.upperNodeIdBound(),std::vector<double>(G.upperNodeIdBound(), 0.0)),
storePreds(storePredecessors) {

}


void DynBetweenness::run() {
/*    if (G.isDirected()) {
        throw std::runtime_error("Invalid argument: G must be undirected.");
    }*/
    count z = G.upperNodeIdBound();
    scoreData.clear();
    scoreData.resize(z);
/*    npaths.clear();
    npaths.resize(z);
    distances.clear();
    distances.resize(z);
    dependencies.clear();
    dependencies.resize(z);
*/
    if (storePreds) {
        predecessors.clear();
        predecessors.resize(G.upperNodeIdBound());
        G.forNodes([&](node v){
            predecessors[v].resize(G.upperNodeIdBound());
        });
    }

    auto computeDependencies = [&](node s) {
        // run SSSP algorithm and keep track of everything
        std::unique_ptr<SSSP> sssp;
        if (G.isWeighted()) {
            sssp.reset(new Dijkstra(G, s, true, true));
        } else {
            sssp.reset(new BFS(G, s, true, true));
        }

        sssp->run();

        G.forNodes([&](node t){
            distances[s][t] = sssp->distance(t);
            npaths[s][t] = sssp->numberOfPaths(t);
            if (storePreds) {
                predecessors[s][t] = sssp->getPredecessors(t);
            }
        });

        // compute dependencies for nodes in order of decreasing distance from s
        std::stack<node> stack = sssp->getStack();
        // set maxDistance to the distance of the furthest vertex
        maxDistance[s] = distances[s][stack.top()];
        while (!stack.empty()) {
            node t = stack.top();
            stack.pop();
            for (node p : sssp->getPredecessors(t)) {
                // workaround for integer overflow in large graphs
                bigfloat tmp = npaths[s][p] / npaths[s][t];
                double weight;
                tmp.ToDouble(weight);

                dependencies[s][p] += weight * (1 + dependencies[s][t]);
            }
            if (t != s) {
                scoreData[t] += dependencies[s][t];
            }
        }
    };

    G.forNodes(computeDependencies);
}


void DynBetweenness::updateUnweighted(GraphEvent e) {
    node u_l, u_h;
    if (G.isDirected()) {
        u_h = e.u;
        u_l = e.v;
    }
    G.forNodes([&] (node s){
        if (!G.isDirected()) {
            if (distances[s][e.u] > distances[s][e.v]){
                u_l = e.u;
                u_h = e.v;
            } else {
                u_l = e.v;
                u_h = e.u;
            }
        }
        edgeweight difference = distances[s][u_l] - distances[s][u_h];
        if (difference > 0) {
            std::vector<bigfloat> new_npaths(G.upperNodeIdBound());
            std::vector<edgeweight> new_dist(G.upperNodeIdBound());
            G.forNodes([&] (node v){
                new_npaths[v] = npaths[s][v];
                new_dist[v] = distances[s][v];
            });
            std::vector<double> new_dep(G.upperNodeIdBound(), 0.0);
            std::vector<int> touched(G.upperNodeIdBound(), 0);
            std::vector<std::queue<node>> l_queues(maxDistance[s]+1);
            std::queue<node> queue_BFS;
            queue_BFS.push(u_l);
            // one-level update
            if (difference == 1) {
                l_queues[distances[s][u_l]].push(u_l);
                if (storePreds)
                    predecessors[s][u_l].push_back(u_h);
                new_npaths[u_l] += npaths[s][u_h];
                touched[u_l] = -1;
                // BFS traversal from u_l
                while (!queue_BFS.empty()) {
                    node v = queue_BFS.front();
                    DEBUG("extracted node ", v);
                    queue_BFS.pop();
                    G.forNeighborsOf(v, [&](node w) {
                        if (new_dist[w] == new_dist[v]+1) {
                            if (touched[w] == 0) {
                                touched[w] = -1;
                                queue_BFS.push(w);
                                l_queues[new_dist[w]].push(w);
                            }
                            new_npaths[w] = new_npaths[w] - npaths[s][v] + new_npaths[v];
                        }
                    });
                }
            } else if (difference > 1) {
                new_dist[u_l] = distances[s][u_h] + 1;
                l_queues[new_dist[u_l]].push(u_l);
                // BFS traversal from u_l
                while (!queue_BFS.empty()) {
                    node v = queue_BFS.front();
                    DEBUG("extracted node ",v);
                    if (storePreds) {
                        predecessors[s][v].clear();
                    }
                    queue_BFS.pop();
                    touched[v] = -1;
                    new_npaths[v] = 0;
                    G.forInNeighborsOf(v, [&] (node w){
                        if (new_dist[w] + 1 == new_dist[v]) {
                            new_npaths[v] += new_npaths[w];
                            if (storePreds) {
                                predecessors[s][v].push_back(w);
                            }
                        }
                    });
                    G.forNeighborsOf(v, [&] (node w){
                        if (new_dist[w] >= new_dist[v] && touched[w] == 0) {
                            if (new_dist[w] > new_dist[v])
                                new_dist[w] = new_dist[v]+1;
                            touched[w] = -1;
                            l_queues[new_dist[w]].push(w);
                            queue_BFS.push(w);
                        }
                    });
                }
            }
            // dependencies accumulation
            DEBUG("Dependency accumulation");
            count level = maxDistance[s];
            while(level > 0) {
                while(!l_queues[level].empty()) {
                    node w = l_queues[level].front();
                    DEBUG("Node ",w);
                    l_queues[level].pop();
                    auto updateDependency = [&](node w, node v) {
                        if (new_dist[v] < new_dist[w]) {
                            if (touched[v] == 0) {
                                touched[v] = 1;
                                new_dep[v] = dependencies[s][v];
                                l_queues[level-1].push(v);
                            }

                            bigfloat tmp = new_npaths[v] / new_npaths[w];
                            double weight;
                            tmp.ToDouble(weight);
                            double new_contrib = weight * (1 + new_dep[w]);
                            new_dep[v] += new_contrib;
                            tmp = npaths[s][v] / npaths[s][w];
                            tmp.ToDouble(weight);
                            double old_contrib = weight * (1 + dependencies[s][w]);
                            if (touched[v] == 1 && (v != u_h or w!= u_l))
                                new_dep[v] -= old_contrib;
                            DEBUG("Parent ", v);
                        }
                    };
                    if (storePreds) {
                        for (node v : predecessors[s][w]) {
                            updateDependency(w, v);
                        }
                    }
                    else {
                        G.forInNeighborsOf(w, [&](node v){
                            updateDependency(w, v);
                        });
                    }
                    if (w != s) {
                        scoreData[w] = scoreData[w] + new_dep[w] - dependencies[s][w];
                    }
                }
                level = level - 1;
            }
            // data structures update
            G.forNodes([&] (node r){
                distances[s][r] = new_dist[r];
                npaths[s][r] = new_npaths[r];
                if (touched[r] != 0)
                    dependencies[s][r] = new_dep[r];
            });
        }
/*

        BFS sssp(G, s, true, true);
        sssp.run();
        G.forNodes([&](node t){
            if(distances[s][t] != sssp.distance(t) || npaths[s][t] != sssp.numberOfPaths(t)) {
                std::cout<<"Source: "<<s<<"  Node "<<t<<std::endl;
                std::cout<<"COMPUTED DISTANCE: "<<distances[s][t]<<"  REAL DISTANCE: "<<sssp.distance(t)<<std::endl;
                std::cout<<"COMPUTED NUMBER OF PATHS: "<<npaths[s][t]<<"  REAL NUMBER OF PATHS: "<<sssp.numberOfPaths(t)<<std::endl;
            }
        });


        // compute dependencies for nodes in order of decreasing distance from s
        std::vector<double> dep(G.upperNodeIdBound());
        std::stack<node> stack = sssp.getStack();
        while (!stack.empty()) {
            node t = stack.top();
            stack.pop();
            for (node p : sssp.getPredecessors(t)) {
                dep[p] += (double(npaths[s][p]) / npaths[s][t])  * (1 + dep[t]);
            }
        }

        G.forNodes([&](node t){
            if((dependencies[s][t]-dep[t]>0.000001 || dependencies[s][t]-dep[t]<-0.000001) && s!=t) {
                std::cout<<"Source: "<<s<<"  Node "<<t<<std::endl;
                std::cout<<"COMPUTED DEPENDENCY: "<<dependencies[s][t]<<"  REAL DEPENDENCY: "<<dep[t]<<std::endl;
            }
        });

*/

    });
}


void DynBetweenness::updateWeighted(GraphEvent e) {
    G.forNodes([&] (node s){
        scoreData[s] = 0.0;
    });
    G.forNodes([&] (node s){
        // update of distances and number of shortest paths
        Aux::PrioQueue<double, node> S(G.upperNodeIdBound());
        G.forNodes([&] (node t){
            auto updatePaths = [&](node u, node v, edgeweight w) {
                if (distances[s][t] == distances[s][u] + w + distances[v][t]){
                    npaths[s][t] = npaths[s][t] + npaths[s][u] * npaths[v][t];
                    if (storePreds) {
                        if (t != v){
                            DEBUG("Adding predecessors from v to t to predecessors from s to t");
                            //predecessors[s][t].insert( predecessors[s][t].end(), predecessors[v][t].begin(), predecessors[v][t].end());
                            for (node pv : predecessors[v][t]) {
                                bool wasAlreadyPred = false;
                                for (node ps : predecessors[s][t]) {
                                    if (ps == pv)
                                        wasAlreadyPred = true;
                                }
                                if (!wasAlreadyPred)
                                    predecessors[s][t].push_back(pv);
                            }
                            DEBUG("Finished adding");
                        }
                        else {
                            DEBUG("Adding u to the predecessors of v");
                            predecessors[s][t].push_back(u);
                            DEBUG("Finished adding");
                        }
                    }
                }
                if (distances[s][t] > distances[s][u] + w + distances[v][t]){
                    distances[s][t] = distances[s][u] + w + distances[v][t];
                    npaths[s][t] = npaths[s][u] * npaths[v][t];
                    if (storePreds) {
                        if (t != v) {
                            DEBUG("Replacing predecessors from v to t with predecessors from s to t");
                            predecessors[s][t] = predecessors[v][t];
                            DEBUG("Finished replacing");
                        }
                        else {
                            DEBUG("Replacing old predecessors of v with u");
                            predecessors[s][t].clear();
                            predecessors[s][t].push_back(u);
                            DEBUG("Finished replacing");
                        }
                    }
                }
            };
            if (s != t) {
                updatePaths(e.u, e.v, e.w);
                if (!G.isDirected())
                    updatePaths(e.v, e.u, e.w);
                S.insert(-1*distances[s][t], t);
                dependencies[s][t] = 0.0;
            }
        });
        // dependency accumulation
        while (S.size() != 0) {
            // extract the node with the maximum distance
            node t = S.extractMin().second;
            if (storePreds) {
                DEBUG("Computing dependencies");
                for (node p : predecessors[s][t]) {
                    bigfloat tmp = npaths[s][p] / npaths[s][t];
                    double weight;
                    tmp.ToDouble(weight);
                    dependencies[s][p] += weight  * (1 + dependencies[s][t]);
                }
            }
            else {
                G.forInNeighborsOf(t, [&] (node p){
                    if (Aux::NumericTools::logically_equal(distances[s][t], distances[s][p] + G.weight(p, t))) {
                        bigfloat tmp = npaths[s][p] / npaths[s][t];
                        double weight;
                        tmp.ToDouble(weight);
                        dependencies[s][p] += weight * (1 + dependencies[s][t]);
                    }
                });
            }
            if (t != s) {
                scoreData[t] += dependencies[s][t];
            }
        }

    });

}

void DynBetweenness::update(GraphEvent e) {
    if (G.isWeighted())
        updateWeighted(e);
    else
        updateUnweighted(e);
}

} /* namespace NetworKit */
