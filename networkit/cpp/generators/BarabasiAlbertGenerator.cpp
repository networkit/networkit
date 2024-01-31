/*
 * BarabasiAlbertGenerator.cpp
 *
 *  Created on: May 28, 2013
 *      Author: forigem, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <atomic>
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_set>
#include <tlx/math/clz.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {
BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k, count nMax, count n0, bool sequential)
    : initGraph(0), k(k), nMax(nMax), sequential(sequential) {
    if (k > nMax)
        throw std::runtime_error("k (number of attachments per node) may not be larger than the "
                                 "number of nodes in the target graph (nMax)");
    if (n0 > nMax)
        throw std::runtime_error("n0 (number of initially connected nodes) may not be larger than "
                                 "the number of nodes in the target graph (nMax)");

    if (n0 < k) {
        if (n0 > 0) {
            WARN("given n0 is smaller than k, setting n0 = k");
        }
        this->n0 = k;
    } else {
        this->n0 = n0;
    }
}

BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k, count nMax, const Graph &initGraph,
                                                 bool sequential)
    : initGraph(initGraph), k(k), nMax(nMax), n0(0), sequential(sequential) {
    if (initGraph.numberOfNodes() != initGraph.upperNodeIdBound())
        throw std::runtime_error("initGraph is expected to have consecutive node ids");
    if (k > nMax)
        throw std::runtime_error("k (number of attachments per node) may not be larger than the "
                                 "number of nodes in the target graph (nMax)");
    if (initGraph.numberOfNodes() > nMax)
        throw std::runtime_error(
            "initialization graph cannot have more nodes than the target graph (nMax)");
}

struct RNG {
    uint64_t seed;
    RNG(uint64_t seed) : seed(seed) {}

    void next(uint64_t &in) { in = (in * seed ^ 58493648962871) * 1432563757536332561L; }

    uint64_t nextInt(uint64_t &in, uint64_t max) {
        const uint64_t shift = tlx::clz(max);
        next(in);
        while (true) {
            next(in);
            uint64_t val = in >> shift;
            if (val < max)
                return val;
        }
    }
};

Graph BarabasiAlbertGenerator::generate() {
    if (!nMax)
        return Graph();

    if (sequential) {
        return generateBatagelj();
    } else {
        return generateParallel();
    }
}

Graph BarabasiAlbertGenerator::generateParallel() {
    const node n = nMax;

    GraphBuilder builder{nMax};
    builder.setAutoCompleteEdges(true);
    Graph G{nMax};

    // Temporarily stored edges of the seed graph in edge list to allow fast random access.
    std::vector<node> seedGraphData;

    auto addSeedGraphEdge = [&](node u, node v) {
        seedGraphData.push_back(u);
        seedGraphData.push_back(v);
        G.addEdge(u, v);
    };

    if (initGraph.numberOfNodes() == 0) {
        seedGraphData.reserve(2 * n0 + 2 * (n - n0) * k);

        // initialize n0 connected nodes
        for (index v = 0; v < n0 - 1; ++v) {
            addSeedGraphEdge(v, v + 1);
        }
        addSeedGraphEdge(0, n0 - 1);
    } else {
        seedGraphData.reserve(2 * initGraph.numberOfEdges()
                              + 2 * (n - initGraph.numberOfNodes()) * k);

        initGraph.forEdges([&](node u, node v) { addSeedGraphEdge(u, v); });
        n0 = initGraph.numberOfNodes();
    }

    // for each of the remaining nodes [n0, n), we draw k random DIFFERENT neighbors
    RNG rng{Aux::Random::integer(0xffffffffffffffffL) | 1};
#pragma omp parallel
    {
        std::vector<index> currentEdges(k);
        int threads = omp_get_num_threads();
        int threadId = omp_get_thread_num();
        for (index v = n0 + threadId; v < n; v += threads) {
            for (index i = 0; i < k; ++i) {
                for (index attempt = 0;; ++attempt) {
                    index u = 2 * (i + k * (v - n0)) + seedGraphData.size();
                    uint64_t seed = u + attempt * 0x52785628791L;
                    while (true) {
                        u = rng.nextInt(seed, u);
                        if (u < seedGraphData.size()) {
                            u = seedGraphData[u];
                            break;
                        } else if ((u - seedGraphData.size()) % 2 != 0) {
                            u = n0 + (u - seedGraphData.size()) / (2 * k);
                            break;
                        }
                        seed = u;
                    }
                    if (v == u)
                        continue;
                    assert(u < v);
                    bool hasEdge = false;
                    for (index j = 0; j < i; j++) {
                        if (currentEdges[j] == u) {
                            hasEdge = true;
                            break;
                        }
                    }
                    if (!hasEdge) {
                        builder.addHalfEdge(u, v);
                        currentEdges[i] = u;
                        break;
                    }
                }
            }
        }
    }
    G = builder.completeGraph();
    return G;
}

Graph BarabasiAlbertGenerator::generateBatagelj() {
    const node n = nMax;

    // Temporarily stored edges in edge list M to allow fast  random access.
    // Degrees are only computed to accelerate graph building later on.
    // TODO: Once we've a fast GraphBuilder remove degree and migrate to GraphBuilder
    std::vector<node> M;
    std::vector<count> degree(n, 0);

    auto addEdge = [&](node u, node v) {
        M.push_back(u);
        M.push_back(v);
        degree[u]++;
        degree[v]++;
    };

    // copy seed graph into M
    if (initGraph.numberOfNodes() == 0) {
        M.reserve(2 * n0 + 2 * (n - n0) * k);

        // initialize n0 connected nodes
        for (index v = 0; v < n0 - 1; ++v) {
            addEdge(v, v + 1);
        }
        addEdge(0, n0 - 1);
    } else {
        M.reserve(2 * initGraph.numberOfEdges() + 2 * (n - initGraph.numberOfNodes()) * k);

        initGraph.forEdges([&](node u, node v) { addEdge(u, v); });
        n0 = initGraph.numberOfNodes();
    }

    // for each of the remaining nodes [n0, n), we draw k random DIFFERENT neighbors
    auto &gen = Aux::Random::getURNG();
    Aux::SignalHandler handler;
    for (index v = n0; v < n; ++v) {
        // If we were to update the range in the next loop, the additionally available nodes
        // would only lead to self-loops are multi-edges.
        std::uniform_int_distribution<size_t> distr{0, M.size() - 1};
        auto firstNeighbor = M.size() + 1;

        for (index i = 0; i < k; ++i) {
            // let's sample a new neighbor and repeat if we're already connected to it
            while (true) {
                const auto randomIndex = distr(gen);
                const auto newNeighbor = M[randomIndex];

                // the last 2*(i-1) positions contain all edges incident to v in to format
                //  Even  Odd   Even  Odd   Even  Odd
                // | v | Neigh | v | Neigh | v | Neigh ...
                // Hence, we need to compare the new neighbor to the previous (i-1) odd positions
                bool alreadyIncident = false;
                for (auto j = firstNeighbor; j < M.size(); j += 2) {
                    assert(M[j] != v); // ensure that we're not off by 1

                    if (M[j] == newNeighbor) {
                        alreadyIncident = true;
                        break;
                    }
                }

                if (!alreadyIncident) {
                    addEdge(v, newNeighbor);
                    break;
                }

                handler.assureRunning();
            }
        }
    }

    Graph G(nMax);

    for (node u = 0; u < n; ++u)
        G.preallocateUndirected(u, degree[u]);

    for (size_t i = 0; i < M.size(); i += 2)
        G.addEdge(M[i], M[i + 1]);

    return G;
}

} /* namespace NetworKit */
