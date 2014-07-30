/*
 * DynBetweenness.h
 *
 *  Created on: 29.07.2014
 *      Author: ebergamini
 */

#include <stack>
#include <queue>
#include <memory>

#include "DynBetweenness.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {

DynBetweenness::DynBetweenness(const Graph& G) : Centrality(G, normalized),
npaths(G.upperNodeIdBound(), std::vector<count>(G.upperNodeIdBound())),
distances(G.upperNodeIdBound(), std::vector<edgeweight>(G.upperNodeIdBound())),
dependencies(G.upperNodeIdBound(),std::vector<double>(G.upperNodeIdBound(), 0.0)) {

}

void DynBetweenness::run() {
    count z = G.upperNodeIdBound();
    scoreData.clear();
    scoreData.resize(z);
    //TODO: resize also npaths, dependencies and distances


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
        });

        // compute dependencies for nodes in order of decreasing distance from s
        std::stack<node> stack = sssp->getStack();
        while (!stack.empty()) {
            node t = stack.top();
            stack.pop();
            for (node p : sssp->getPredecessors(t)) {
                dependencies[s][p] += (double(npaths[s][p]) / npaths[s][t])  * (1 + dependencies[s][t]);
            }
            if (t != s) {
                scoreData[t] += dependencies[s][t];
            }
        }
    };

    G.forNodes(computeDependencies);
}

void DynBetweenness::update(const GraphEvent e) {
    G.forNodes([&] (node s){
        node u_l, u_h;
        if (distances[s][e.u] > distances[s][e.v]){
            u_l = e.u;
            u_h = e.v;
        } else {
            u_l = e.v;
            u_h = e.u;
        }
        edgeweight difference = distances[s][u_l] - distances[s][u_h];
        if (difference > 0) {
            std::vector<count> new_npaths(G.upperNodeIdBound());
            std::vector<edgeweight> new_dist(G.upperNodeIdBound());
            G.forNodes([&] (node v){
                new_npaths[v] = npaths[s][v];
                new_dist[v] = distances[s][v];
            });
            std::vector<double> new_dep(G.upperNodeIdBound(), 0.0);
            std::vector<int> touched(G.upperNodeIdBound(), 0);
            std::vector<std::queue<node>> l_queues(G.upperNodeIdBound());
            std::queue<node> queue_BFS;
            queue_BFS.push(u_l);
            // one-level update
            if (difference == 1) {
                l_queues[distances[s][u_l]].push(u_l);
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
                // dependencies accumulation
                count level = G.upperNodeIdBound() - 1;
                while (level > 0) {
                    while (!l_queues[level].empty()){
                        node w = l_queues[level].front();
                        l_queues[level].pop();
                        G.forNeighborsOf(w, [&](node v){
                            // v is a predecessor of w
                            if (new_dist[v] < new_dist[w]) {
                                if (touched[v] == 0) {
                                    touched[v] = 1;
                                    new_dep[v] = dependencies[s][v];
                                    l_queues[level - 1].push(v);
                                }
                                double new_contrib = double(new_npaths[v])/new_npaths[w]*(1+new_dep[w]);
                                new_dep[v] += new_contrib;
                                double old_contrib = double(npaths[s][v])/npaths[s][w]*(1+dependencies[s][w]);
                                if (touched[v] == 1 && (v != u_h or w!= u_l))
                                    new_dep[v] -= old_contrib;
                            }
                        });
                        if (w != s) {
                            scoreData[w] = scoreData[w] + new_dep[w] - dependencies[s][w];
                        }
                    }
                    level = level - 1;
                }
            } else if (difference > 1) {
                new_dist[u_l] = distances[s][u_h] + 1;
                l_queues[new_dist[u_l]].push(u_l);
                // BFS traversal from u_l
                while (!queue_BFS.empty()) {
                    node v = queue_BFS.front();
                    DEBUG("extracted node ",v);
                    queue_BFS.pop();
                    touched[v] = -1;
                    new_npaths[v] = 0;
                    G.forNeighborsOf(v, [&] (node w){
                            if (new_dist[w] + 1 == new_dist[v])
                                new_npaths[v] += new_npaths[w];
                            if (new_dist[w] >= new_dist[v] && touched[w] == 0) {
                                if (new_dist[w] > new_dist[v])
                                    new_dist[w] = new_dist[v]+1;
                                touched[w] = -1;
                                l_queues[new_dist[w]].push(w);
                                queue_BFS.push(w);
                            }
                    });
                }
                // dependencies accumulation
                DEBUG("Dependency accumulation");
                count level = G.upperNodeIdBound() - 1;
                while(level > 0) {
                    while(!l_queues[level].empty()) {
                        node w = l_queues[level].front();
                        DEBUG("Node ",w);
                        l_queues[level].pop();
                        G.forNeighborsOf(w, [&](node v){
                            if (new_dist[v] < new_dist[w]) {
                                if (touched[v] == 0) {
                                    touched[v] = 1;
                                    new_dep[v] = dependencies[s][v];
                                    l_queues[level-1].push(v);
                                }
                                double new_contrib = double(new_npaths[v])/new_npaths[w]*(1+new_dep[w]);
                                new_dep[v] += new_contrib;
                                double old_contrib = double(npaths[s][v])/npaths[s][w]*(1+dependencies[s][w]);
                                if (touched[v] == 1 && (v != u_h or w!= u_l))
                                    new_dep[v] -= old_contrib;
                                DEBUG("Parent ", v);
                            }
                        });
                        if (w != s) {
                            scoreData[w] = scoreData[w] + new_dep[w] - dependencies[s][w];
                        }
                    }
                    level = level - 1;
                }
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
            if((dependencies[s][t]-dep[t]>0.0001 || dependencies[s][t]-dep[t]<-0.0001) && s!=t) {
                std::cout<<"Source: "<<s<<"  Node "<<t<<std::endl;
                std::cout<<"COMPUTED DEPENDENCY: "<<dependencies[s][t]<<"  REAL DEPENDENCY: "<<dep[t]<<std::endl;
            }
        });

        */

    });

}

} /* namespace NetworKit */
