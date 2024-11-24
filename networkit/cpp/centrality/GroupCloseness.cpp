/*
 * GroupCloseness.cpp
 *
 *  Created on: 03.10.2016
 *      Author: elisabetta bergamini
 */

#include <atomic>
#include <memory>
#include <omp.h>
#include <queue>
#include <ranges>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/centrality/TopCloseness.hpp>
namespace NetworKit {

GroupCloseness::GroupCloseness(const Graph &G, count k, count H) : G(&G), k(k), H(H) {
    d1Global.resize(omp_get_max_threads(), std::vector<count>(G.upperNodeIdBound()));
}

edgeweight GroupCloseness::computeImprovement(node u, count h) {
    // computes the marginal gain due to adding u to S
    auto &d1 = d1Global[omp_get_thread_num()];
    std::ranges::copy(d, d1.begin());

    d1[u] = 0;
    count improvement = d[u]; // old distance of u
    std::queue<node> Q;
    Q.push(u);

    index level = 0;
    while (!Q.empty() && (h == 0 || level <= h)) {
        node v = Q.front();
        Q.pop();
        level = d1[v];
        G->forNeighborsOf(v, [&](node w) {
            if (d1[w] > d1[v] + 1) {
                d1[w] = d1[v] + 1;
                improvement += d[w] - d1[w];
                Q.push(w);
            }
        });
    }
    return improvement;
}

void GroupCloseness::updateDistances(node u) {
    count improvement = d[u]; // old distance of u
    d[u] = 0;                 // now u is in the group
    std::queue<node> Q;
    Q.push(u);

    do {
        node v = Q.front();
        Q.pop();
        G->forNeighborsOf(v, [&](node w) {
            if (d[w] > d[v] + 1) {
                improvement += d[w] - (d[v] + 1);
                d[w] = d[v] + 1;
                Q.push(w);
            }
        });
    } while (!Q.empty());
}

void GroupCloseness::run() {
    const count n = G->upperNodeIdBound();
    node top;

    omp_lock_t lock;
    omp_init_lock(&lock);

    if (H == 0) {
        TopCloseness topcc(*G, 1, true, false);
        topcc.run();
        top = topcc.topkNodesList()[0];
    } else
        top = *std::ranges::max_element(
            G->nodeRange().begin(), G->nodeRange().end(),
            [&G = G](node u, node v) { return G->degree(u) > G->degree(v); });

    // first, we store the distances between each node and the top node
    d.clear();
    d.resize(n);

    Traversal::BFSfrom(*G, top, [&d = d](node u, count distance) { d[u] = distance; });

    count sumD = G->parallelSumForNodes([&](node v) { return d[v]; });

    // init S
    S.clear();
    S.resize(k, 0);
    S[0] = top;

    std::vector<int64_t> prevBound(n, 0);
    count currentImpr = sumD + 1; // TODO change

    std::vector<count> S2(n, sumD);
    S2[top] = 0;

    // loop to find k group members
    for (index i = 1; i < k; i++) {
        DEBUG("k = ", i);
        Aux::BucketPQ Q(n, -currentImpr - 1, 0);
        G->forNodes([&](node v) {
            if (d[v] > 0)
                Q.insert(-prevBound[v], v);
        });
        currentImpr = 0;
        node maxNode = none;

        std::atomic<bool> toInterrupt{false};
#pragma omp parallel
        {
            while (!toInterrupt.load(std::memory_order_relaxed)) {
                omp_set_lock(&lock);
                if (Q.empty()) {
                    omp_unset_lock(&lock);
                    toInterrupt.store(true, std::memory_order_relaxed);
                    break;
                }

                const node v = Q.extractMin().second;

                omp_unset_lock(&lock);
                INFO("Extracted node ", v, " with prio ", prevBound[v]);
                if (i > 1 && prevBound[v] <= static_cast<int64_t>(currentImpr)) {
                    INFO("Interrupting! currentImpr = ", currentImpr,
                         ", previous bound = ", prevBound[v]);
                    toInterrupt.store(true, std::memory_order_relaxed);
                    break;
                }
                if (i == 1 || prevBound[v] > static_cast<int64_t>(currentImpr)) {
                    count imp = computeImprovement(v, H);
                    omp_set_lock(&lock);
                    if (imp > currentImpr) {
                        currentImpr = imp;
                        maxNode = v;
                    }
                    omp_unset_lock(&lock);
                    INFO("New bound for v = ", imp);
                    prevBound[v] = imp; // TODO use later
                } else {
                    prevBound[v] = 0;
                }
            }
        }
        S[i] = maxNode;

        updateDistances(S[i]);
    }

    hasRun = true;
}

double GroupCloseness::computeFarness(const std::vector<node> &S, count H) const {
    // we run a BFS from S up to distance H (if H > 0) and sum the distances
    double farness = 0;
    count k = S.size();
    std::vector<double> d1(G->upperNodeIdBound(), 0);
    std::vector<bool> visited(G->upperNodeIdBound(), false);
    std::queue<node> Q;

    for (node i = 0; i < k; i++) {
        Q.push(S[i]);
        visited[S[i]] = true;
    }

    while (!Q.empty()) {
        node u = Q.front();
        Q.pop();
        if (H > 0 && d1[u] > H) {
            break;
        }
        farness += d1[u];
        G->forNeighborsOf(u, [&](node w) {
            if (!visited[w]) {
                visited[w] = true;
                d1[w] = d1[u] + 1;
                Q.push(w);
            }
        });
    }
    return farness;
}

} /* namespace NetworKit */
