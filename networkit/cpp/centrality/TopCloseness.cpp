/*
 * TopCloseness.cpp
 *
 *  Created on: 03.42658.2014
 *      Author: nemes
 */

#include <stack>
#include <queue>
#include <memory>
#include <omp.h>

#include "TopCloseness.h"
#include "../components/ConnectedComponents.h"
#include "../components/StronglyConnectedComponents.h"
#include "../auxiliary/PrioQueueForInts.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {

TopCloseness::TopCloseness(const Graph& G, count k, bool first_heu, bool sec_heu) : G(G), k(k), first_heu(first_heu), sec_heu(sec_heu) {
}


void TopCloseness::init() {
    n = G.upperNodeIdBound();
    topk.clear();
    topk.resize(k);
    topkScores.clear();
    topkScores.resize(k);
    DEBUG("Number of nodes: ", n);
    DEBUG("k = ", k);
    farness.clear();
    farness.resize(n, 0);
    computeReachable();
    DEBUG("Done INIT");
}

void TopCloseness::computeReachable() {
    if (G.isDirected()) {
        computeReachableNodesDir();
    } else {
        computeReachableNodesUndir();
    }
}

void TopCloseness::computeReachableNodesDir() {
    reachL = std::vector<count>(n);
    reachU = std::vector<count>(n);
    component = std::vector<count>(n);
    DEBUG("Before running SCCs");
    StronglyConnectedComponents sccs(G);
    sccs.run();

    count N = sccs.numberOfComponents();
    DEBUG("Number of components: ", N);
    std::vector<count> reachL_scc(N, 0);
    std::vector<count> reachU_scc(N, 0);
    std::vector<count> reachU_without_max_scc(N, 0);
    std::vector<bool> reach_from_max_scc(N, false);
    std::vector<bool> reaches_max_scc(N, false);
    std::vector<std::vector<count>> sccs_vec(N, std::vector<count>());
    Graph sccGraph(N, false, true);
    std::vector<bool> found(N, false);
    count maxSizeCC = 0;

    // We compute the vector sccs_vec, where each component contains the list of its nodes
    for (count v = 0; v < n; v++) {
        component[v] = sccs.componentOfNode(v);
        sccs_vec[sccs.componentOfNode(v)-1].push_back(v);
    }

    // We compute the SCC graph and store it in sccGraph
    for (count V = 0; V < N; V++) {
        for (count v:sccs_vec[V]) {
            G.forNeighborsOf(v, [&](node w){
                count W = sccs.componentOfNode(w) - 1;

                if (W != V && !found[W]) {
                    found[W] = true;
                    sccGraph.addEdge(V, W);
                }
            });
        }
        sccGraph.forNeighborsOf(V, [&](node W) {
            found[W] = false;
        });
        if (sccGraph.degreeOut(V) > sccGraph.degreeOut(maxSizeCC)) {
            maxSizeCC = V;
        }
        // ELISABETTA: maybe the code can be made simpler by running G.forEdges to scan all the edges. Would it be better to have a Graph object to store the SCC graph?
    }     // MICHELE: I have used a graph instead of scc_adjlist. About G.forEdges, I think it is
    // a bit more complicated: I have to scan nodes, otherwise I do not know how to avoid multiple edges.
    // This scan is made using variable "found". Do you have better ideas? Note that this is linear in the graph
    // size.

    //BFS from the biggest SCC.
    std::queue<count> Q;
    Q.push(maxSizeCC);
    reach_from_max_scc[maxSizeCC] = true;
    while (!Q.empty()) {
        count V = Q.front();
        Q.pop();
        reachL_scc[maxSizeCC] += sccs_vec[V].size();
        sccGraph.forNeighborsOf(V, [&](node W) {
            if (!reach_from_max_scc[W]) {
                reach_from_max_scc[W] = true;
                Q.push(W);
            }
        });
    }
    reachU_scc[maxSizeCC] = reachL_scc[maxSizeCC];
    reaches_max_scc[maxSizeCC] = true;

    //so far only the largest SCC has reach_U and reach_L > 0

    // Dynamic programming to compute number of reachable vertices
    for (count V = 0; V < N; V++) {
        if (V == maxSizeCC) {
            continue;
        }
        sccGraph.forNeighborsOf(V, [&](node W) {
            reachL_scc[V] = std::max(reachL_scc[V], reachL_scc[W]);
            if (!reach_from_max_scc[W]) {
                reachU_without_max_scc[V] += reachU_without_max_scc[W];
            }
            reachU_scc[V] += reachU_scc[W];
            reachU_scc[V] = std::min(reachU_scc[V], n);
            reaches_max_scc[V] = reaches_max_scc[V] || reaches_max_scc[W];
        });

        if (reaches_max_scc[V]) {
            reachU_scc[V] = reachU_without_max_scc[V] + reachU_scc[V];
        }
        reachL_scc[V] += sccs_vec[V].size();
        reachU_scc[V] += sccs_vec[V].size();
        reachU_scc[V] = std::min(reachU_scc[V], n);
    }

    for (count v = 0; v < n; v++) {
        reachL[v] = reachL_scc[sccs.componentOfNode(v)-1];
        reachU[v] = reachU_scc[sccs.componentOfNode(v)-1];
        if (false) { // MICHELE: used to check if the bounds are correct
            count r = 0;
            G.BFSfrom(v,[&](node w, count dist) {
                r++;
            });

            if (reachL[v] > r || reachU[v] < r) {
                DEBUG("BIG MISTAKE! ", reachL[v], " ", r, " ", reachU[v]);
                while(true) {}
            }
        }
    }
}

void TopCloseness::computeReachableNodesUndir() {
    reachL = std::vector<count>(n);

    ConnectedComponents comps(G);
    comps.run();
    std::map<index, count> sizes = comps.getComponentSizes();
    G.forNodes([&](node v){
        index cv = comps.componentOfNode(v);
        reachL[v] = sizes[cv];
    });
    reachU = reachL;
}

void TopCloseness::computelBound1(std::vector<double> &S) {
    std::vector<count> neighbors(n, 0);
    std::vector<count> N(n, 0);
    std::vector<count> neighbors_new(n, 0);
    std::vector<count> neighbors_old(n, 0);
    std::vector<double> sumDist(n, 0);
    std::vector<bool> finished(n, false);

    count n_finished = 0;


    G.forNodes([&](node u){
        S[u] = std::numeric_limits<double>::max();
        if (G.degreeOut(u) == 0) {
            finished[u] = true;
            n_finished++;
        }
        neighbors[u] = G.degreeOut(u);
        sumDist[u] = neighbors[u];
        N[u] = neighbors[u] + 1; // we also count the node itself in the number of visited nodes
    });
    count level = 2;
    DEBUG("computing first lbound");


    while(n_finished < n) {
        DEBUG("First bound. Finished: ", n_finished, " of ", n, ".");
        G.forNodes([&](node u){
            if (!finished[u]) {
                n_op += G.degreeOut(u);
                neighbors_new[u] = 0;
                G.forNeighborsOf(u, [&](node v){
                    neighbors_new[u] += neighbors[v];
                });
                if (!G.isDirected()) {
                    if (level == 2) {
                        neighbors_new[u] -= G.degreeOut(u);
                    } else {
                        if (neighbors_new[u] < (G.degreeOut(u) - 1) * neighbors_old[u]) {
                            DEBUG("BIG MISTAKE");
                            while(true) {}
                        }
                        neighbors_new[u] -= (G.degreeOut(u) - 1) * neighbors_old[u];
                    }
                }

                count n_old = N[u];
                N[u] += neighbors_new[u];
                sumDist[u] += level*neighbors_new[u];

                if (N[u] >= reachL[u]) {
                    if (n_old < reachL[u]) {
                        // We have to consider the case in which the number of reachable vertices is reachL.
                        S[u] = (sumDist[u] - level*(N[u]-reachL[u])) * (n-1) / (reachL[u]-1) / (reachL[u]-1);
                    }
                    if (neighbors_new[u] == 0) {
                        reachU[u] = N[u];
                    }
                    if (N[u] >= reachU[u]) {
                        // We have to consider the case in which the number of reachable vertices is reachU.
                        S[u] = std::min(S[u], (sumDist[u] - level*(N[u]-reachU[u])) * (n-1) / (reachU[u]-1) / (reachU[u]-1));
                        finished[u] = true;
                        n_finished++;

                        if (N[u] < reachL[u] && neighbors_new[u] == 0) {
                            DEBUG("BIG MISTAKE!!!", reachL[u]);
                        }

                    } else { // reachL < N < reachU
                        // We have to consider the case in which the number of reachable is N[u].
                        S[u] = std::min(S[u], sumDist[u] * (n-1) / (N[u]-1) / (N[u]-1));
                    }
                }
            }
        });
        G.forNodes([&](node u){
            // We update neighbors.
            neighbors_old[u] = neighbors[u];
            neighbors[u] = neighbors_new[u];
        });
        level ++;
    }
    DEBUG("Visited edges (first lbound): ", n_op);
}

void TopCloseness::BFSbound(node x, std::vector<double> &S2, count *visEdges) {
    count r = 0;
    std::vector<std::vector<node>> levels(n);
    // nodesPerLev[i] contains the number of nodes in level i
    std::vector<count> nodesPerLev(n, 0);
    // sumLevs[i] contains the sum of the nodes in levels j <= i
    std::vector<count> sumLevs(n, 0);
    count nLevs = 0;
    levels[nLevs].clear();
    double sum_dist = 0;
    G.BFSfrom(x,[&](node u, count dist) {
        sum_dist += dist;
        r++;
        if (dist > nLevs) {
            sumLevs[nLevs] += nodesPerLev[nLevs];
            sumLevs[nLevs+1] = sumLevs[nLevs];
            nLevs ++;
            levels[nLevs].clear();
        }
        levels[nLevs].push_back(u);
        nodesPerLev[nLevs] ++;
    });
    sumLevs[nLevs] += nodesPerLev[nLevs];
    if (G.isDirected()) {
        (*visEdges) += G.numberOfEdges();
    } else {
        (*visEdges) += 2 * G.numberOfEdges();
    }
    S2[x] = sum_dist * (n - 1.0) / (r - 1.0) / (r - 1.0);
    // we compute the bound for the first level
    count closeNodes = 0, farNodes = 0;
    for (count j = 0; j <= nLevs; j++) {
        if (abs((long long)j-1LL)<=1)  {
            closeNodes += nodesPerLev[j];
        } else {
            farNodes += nodesPerLev[j]*abs(1LL-(long long)j);
        }
    }

    edgeweight level_bound = 2.0*(closeNodes) + (double)farNodes;
    for (count j = 0; j < levels[1].size(); j ++) {
        node w = levels[1][j];
        // we subtract 2 not to count the node itself
        double bound = (level_bound - 2 - G.degree(w)) * (n-1.0) / (reachU[w]-1.0) / (reachU[w]-1.0);
        if (bound > S2[w] && (!G.isDirected() || component[w] == component[x])) {
            S2[w] = bound;
        }
    }
    //DEBUG("level_bound = ", level_bound);
    // now we compute it for the other levels
    for (count i = 2; i <= nLevs; i++) {
        if (!G.isDirected() && i > 2) {
            level_bound +=sumLevs[i-3];
        }
        if (i < nLevs) {
            level_bound -= (sumLevs[nLevs] - sumLevs[i+1]);
        }
        for (count j = 0; j < levels[i].size(); j ++) {
            node w = levels[i][j];
            double bound = (level_bound - 2 - G.degree(w)) * (n-1.0) / (reachU[w]-1.0) / (reachU[w]-1.0);
            if (bound > S2[w] && (!G.isDirected() || component[w] == component[x])) {
                // TODO MICHELE: as before.
                S2[w] = bound;
            }
        }
    }
}



double TopCloseness::BFScut(node v, double x, bool *visited, count *distances, node *pred, count *visEdges) {
    count d = 0, f = 0, nd = 1;
    double rL = reachL[v], rU = reachU[v];
    std::queue<node> Q1;
    std::queue<node> to_reset;
    count sum_dist = 0;
    double ftildeL=0, ftildeU=0, gamma=G.degreeOut(v);
    double farnessV = 0;

    // MICHELE: variable visited is not local, otherwise the allocation would be too expensive.
    visited[v] = true;
    distances[v] = 0;
    Q1.push(v);
    to_reset.push(v);

    while (!Q1.empty()) {
        node u = Q1.front();
        Q1.pop();

        sum_dist += distances[u];
        if (distances[u] > d) { // Need to update bounds!
            d++;
            ftildeL = (f + (d+2)*(rL-nd) - gamma) * (n-1) / (rL-1.0) / (rL-1.0);
            ftildeU = (f + (d+2)*(rU-nd) - gamma) * (n-1) / (rU-1.0) / (rU-1.0);
            if (std::min(ftildeL, ftildeU) >= x) {
                farnessV = std::min(ftildeL, ftildeU);
                break;
            }
            gamma = 0;
        }
        bool cont = true;
        G.forNeighborsOf(u, [&](node w){
            if (cont) {
                (*visEdges)++;
                if (!visited[w]) {
                    distances[w] = distances[u] + 1;
                    Q1.push(w);
                    to_reset.push(w);
                    visited[w] = true;
                    f += distances[w];
                    if (!G.isDirected())
                        gamma += (G.degree(w) - 1); // notice: only because it's undirected
                    else
                        gamma += G.degreeOut(w);
                    nd ++;
                    pred[w] = u;
                } else {
                    if (G.isDirected() || pred[u]!= w) {
                        ftildeL += (n-1) / (rL-1.0) / (rL-1.0);
                        ftildeU += (n-1) / (rU-1.0) / (rU-1.0);
                        if (std::min(ftildeL, ftildeU) >= x) {
                            cont = false;
                        }
                    }
                }
            }
        });
        if (std::min(ftildeL, ftildeU) >= x) {
            farnessV = std::min(ftildeL, ftildeU);
            break;
        }
    }
    while(!to_reset.empty()) { // MICHELE: need to reset variable visited. Variables pred and distances
        // do not need to be updated.
        node u = to_reset.front();
        to_reset.pop();
        visited[u] = false;
    }
    if (farnessV < x) {
        farnessV = sum_dist * (n-1) / (nd-1.0) / (nd-1.0);
    }
    return farnessV;
}



void TopCloseness::run() {
    init();
    Aux::PrioQueue<double, node> top(n); // like in Crescenzi
    std::vector<bool> toAnalyze(n, true);
    omp_lock_t lock;
    omp_init_lock(&lock);

    std::vector<double> S(n);
    // first lower bound on s
    if (first_heu) {
        DEBUG("Computing Neighborhood-based lower bound");
        computelBound1(S);
    }

    DEBUG("Initializing queue");
    G.forNodes([&](node u){
        if (G.degreeOut(u) == 0) {
            farness[u] = std::numeric_limits<double>::max();
        } else if (first_heu) {
            farness[u] = S[u];
        } else {
            farness[u] = -((double) G.degreeOut(u));
        }
    });
    Aux::PrioQueue<double, node> Q(farness);
    DEBUG("Done filling the queue");

#pragma omp parallel // Shared variables:
    // cc: synchronized write, read leads to a positive race condition;
    // Q: fully synchronized;
    // top: fully synchronized;
    // toAnalyze: fully synchronized;
    // visEdges: one variable for each thread, summed at the end;
    {
        double kth = std::numeric_limits<double>::max(); //like in Crescenzi
        bool *visited = NULL;
        count *distances = NULL;
        node *pred = NULL;
        count visEdges = 0;
        #if LOG_LEVEL >= LOG_LEVEL_DEBUG
        count iters = 0;
        #endif
        double farnessS;

        if (omp_get_thread_num() == 0) {
            DEBUG("Number of threads: ", omp_get_num_threads());
        }

        if (!sec_heu) {
            visited = (bool*) calloc(n, sizeof(bool));
            distances = (count*) malloc(n * sizeof(count));
            pred = (node*) malloc(n * sizeof(node));
        }

        while(Q.size() != 0) {
            DEBUG("To be analyzed: ", Q.size());
            omp_set_lock(&lock);
            if (Q.size() == 0) { // The size of Q might have changed.
                omp_unset_lock(&lock);
                break;
            }
            std::pair<double, node> p = Q.extractMin(); // Access to Q must be synchronized
            node s = p.second; //TODO change!
            toAnalyze[s] = false;
            omp_unset_lock(&lock);

            if (G.degreeOut(s) == 0 || farness[s] >= kth) {
                break;
            }
            DEBUG("Iteration ", ++iters, " of thread ", omp_get_thread_num());

            DEBUG("    Extracted node ", s, " with priority ", p.first, ".");
            if (G.degreeOut(s) == 0) {

                omp_set_lock(&lock);
                toAnalyze[s] = false;
                farness[s] = std::numeric_limits<double>::max();
                omp_unset_lock(&lock);

            } else if (sec_heu) {
                // MICHELE: we use BFSbound to bound the centrality of all nodes.
                DEBUG("    Running BFSbound.");
                BFSbound(s, S, &visEdges);
                omp_set_lock(&lock);
                farness[s] = S[s];
                omp_unset_lock(&lock);
                count imp = 0;
                for (count v = 0; v < n; v++) {
                    if (farness[v] < S[v] && toAnalyze[v]) { // This part must be syncrhonized.
                        omp_set_lock(&lock);
                        if (farness[v] < S[v] && toAnalyze[v]) { // Have to check again, because the variables might have changed
                            imp++;
                            farness[v] = S[v];
                            Q.remove(v);
                            Q.insert(farness[v], v);
                        }
                        omp_unset_lock(&lock);
                    }
                }
                DEBUG("    We have improved ", imp, " bounds.");
            } else {
                // MICHELE: we use BFScut to bound the centrality of s.
                DEBUG("    Running BFScut with x=", kth, " (degree:", G.degreeOut(s), ").");
                farnessS = BFScut(s, kth, visited, distances, pred, &visEdges);
                DEBUG("    Visited edges: ", visEdges, ".");
                omp_set_lock(&lock);
                farness[s] = farnessS;
                omp_unset_lock(&lock);
            }

            // If necessary, we update kth.
            if (farness[s] < kth) {
                omp_set_lock(&lock);
                DEBUG("    The closeness of s is ", 1.0 / farness[s], ".");
                top.insert(-farness[s], s);
                if (top.size() > k) {
                    top.extractMin();
                }
                omp_unset_lock(&lock);
            } else {
                DEBUG("    Not in the top-k.");
            }

            // We load the new value of kth.

            if (top.size() == k) {
                omp_set_lock(&lock);
                std::pair<double, node> elem = top.extractMin();
                kth = -elem.first;
                top.insert(elem.first, elem.second);
                omp_unset_lock(&lock);
            }
        }
        DEBUG("Number of iterations of thread ", omp_get_thread_num(), ": ", iters, " out of ", n);
        if (!sec_heu) {
            free(visited);
            free(distances);
            free(pred);
        }
        omp_set_lock(&lock);
        this->visEdges += visEdges;
        omp_unset_lock(&lock);

    }

    hasRun = true;
    for (int i = top.size() - 1; i >= 0; i--) {
        std::pair<double, node> elem = top.extractMin();
        topk[i] = elem.second;
        topkScores[i] = 1.0 / -elem.first;
    }
    for(count j = 0; j < k; j ++) {
        DEBUG(j+1,"-th node with max closeness: ", topk[j], ", its closeness: ", topkScores[j]);
    }
}

} /* namespace NetworKit */
