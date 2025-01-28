/*
 * GraphBuilder.cpp
 *
 *  Created on: 15.07.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <omp.h>
#include <stdexcept>

#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

GraphBuilder::GraphBuilder(count n, bool weighted, bool directed, bool autoCompleteEdges)
    : n(n), selfloops(0), weighted(weighted), directed(directed),
      autoCompleteEdges(autoCompleteEdges) {
    // pre-allocating the adjacency vectors
    const int max_threads = omp_get_max_threads();
    outEdgesPerThread.resize(max_threads);
    inEdgesPerThread.resize(max_threads);
    for (int i = 0; i < max_threads; ++i) {
        outEdgesPerThread[i].resize(max_threads);
        inEdgesPerThread[i].resize(max_threads);
    }
    if (weighted) {
        outEdgeWeightsPerThread.resize(max_threads);
        inEdgeWeightsPerThread.resize(max_threads);
        for (int i = 0; i < max_threads; ++i) {
            outEdgeWeightsPerThread[i].resize(max_threads);
            inEdgeWeightsPerThread[i].resize(max_threads);
        }
    }
}

void GraphBuilder::reset(count n) {
    this->n = n;
    selfloops = 0;

    const int max_threads = omp_get_max_threads();
    outEdgesPerThread.assign(max_threads, std::vector<std::vector<HalfEdge>>{});
    outEdgeWeightsPerThread.assign(isWeighted() ? max_threads : 0,
                                   std::vector<std::vector<edgeweight>>{}),
        inEdgesPerThread.assign(isDirected() ? max_threads : 0,
                                std::vector<std::vector<HalfEdge>>{}),
        inEdgeWeightsPerThread.assign((isDirected() && isWeighted()) ? max_threads : 0,
                                      std::vector<std::vector<edgeweight>>{});
    outEdgesPerThread.resize(max_threads);
    for (int i = 0; i < max_threads; ++i) {
        outEdgesPerThread[i].resize(max_threads);
    }
    if (weighted) {
        outEdgeWeightsPerThread.resize(max_threads);
        for (int i = 0; i < max_threads; ++i) {
            outEdgeWeightsPerThread[i].resize(max_threads);
        }
    }
    inEdgesPerThread.resize(max_threads);
    if (weighted)
        inEdgeWeightsPerThread.resize(max_threads);
    for (int i = 0; i < max_threads; ++i) {
        inEdgesPerThread[i].resize(max_threads);
        if (weighted)
            inEdgeWeightsPerThread[i].resize(max_threads);
    }
}

index GraphBuilder::indexInOutEdgeArrayPerThread(node u, node v) const {
    int max_threads = omp_get_max_threads();
    int thread_num = omp_get_thread_num();
    auto &edges = outEdgesPerThread[thread_num][u % max_threads];
    index result = none;
    for (index i = 0; i < edges.size(); ++i) {
        if (edges[i].source == u && edges[i].destination == v) {
            result = i;
            break;
        }
    }
    return result;
}

index GraphBuilder::indexInInEdgeArrayPerThread(node u, node v) const {
    int max_threads = omp_get_max_threads();
    int thread_num = omp_get_thread_num();
    auto &edges = inEdgesPerThread[thread_num][u % max_threads];
    index result = none;
    for (index i = 0; i < edges.size(); ++i) {
        if (edges[i].source == u && edges[i].destination == v) {
            result = i;
            break;
        }
    }
    return result;
}

node GraphBuilder::addNode() {
    return n++;
}

void GraphBuilder::addHalfEdge(index a, index b, edgeweight ew) {
    if (autoCompleteEdges) {
        if (directed) {
            addHalfOutEdge(a, b, ew);
            addHalfInEdge(b, a);
        } else {
            addHalfOutEdge(a, b, ew);
            if (a != b)
                addHalfOutEdge(b, a, ew);
        }
    } else {
        addHalfOutEdge(a, b, ew);
    }
    if (a == b) {
#pragma omp atomic
        selfloops++;
    }
}

void GraphBuilder::addHalfOutEdge(node u, node v, edgeweight ew) {
    outEdgesPerThread[omp_get_thread_num()][u % omp_get_max_threads()].emplace_back(u, v);
    if (weighted) {
        outEdgeWeightsPerThread[omp_get_thread_num()][u % omp_get_max_threads()].emplace_back(ew);
    }
}

void GraphBuilder::addHalfInEdge(node u, node v, edgeweight ew) {
    inEdgesPerThread[omp_get_thread_num()][u % omp_get_max_threads()].emplace_back(u, v);
    if (weighted) {
        inEdgeWeightsPerThread[omp_get_thread_num()][u % omp_get_max_threads()].emplace_back(ew);
    }
}

void GraphBuilder::swapNeighborhood(node u, std::vector<node> &neighbours,
                                    std::vector<edgeweight> &weights, bool selfloop) {
    if (weighted)
        assert(neighbours.size() == weights.size());

    const index max_threads = omp_get_max_threads();

    std::vector<HalfEdge> new_edges(neighbours.size());
    for (index i = 0; i < neighbours.size(); ++i) {
        new_edges[i] = HalfEdge(u, neighbours[i]);
    }

    const size_t size_per_thread = static_cast<int>(new_edges.size() / max_threads);
    const size_t remainder = new_edges.size() % max_threads;
    std::vector<std::vector<HalfEdge>> new_edges_per_thread(max_threads);

    // split new neighbours into almost equally distributed parts
    int cur_index = 0;
    for (index cur_thread = 0; cur_thread < max_threads; ++cur_thread) {
        const size_t count = size_per_thread + (cur_thread < remainder ? 1 : 0);
        new_edges_per_thread[cur_thread].assign(new_edges.begin() + cur_index,
                                                new_edges.begin() + cur_index + count);
        cur_index += count;
    }
    // distribute new edges across all threads
    for (index thread = 0; thread < max_threads; ++thread) {
        auto current_edges = &outEdgesPerThread[thread][u % max_threads];
        current_edges->swap(new_edges_per_thread[thread]);
    }

    if (weighted) {
        std::vector<std::vector<edgeweight>> new_weights_per_thread(max_threads);

        // split new weights into almost equally distributed parts
        int cur_index_weighted = 0;
        for (index cur_thread = 0; cur_thread < max_threads; ++cur_thread) {
            const size_t count = size_per_thread + (cur_thread < remainder ? 1 : 0);
            new_weights_per_thread[cur_thread].assign(weights.begin() + cur_index_weighted,
                                                      weights.begin() + cur_index_weighted + count);
            cur_index_weighted += count;
        }
        // distribute new weights across all threads
        for (index thread = 0; thread < max_threads; ++thread) {
            auto current_weights = &outEdgeWeightsPerThread[thread][u % max_threads];
            current_weights->swap(new_weights_per_thread[thread]);
        }
    }

    if (selfloop) {
#pragma omp atomic
        selfloops++;
    }
}
void GraphBuilder::setOutWeight(node u, node v, edgeweight ew) {
    assert(isWeighted());
    const index vi = indexInOutEdgeArrayPerThread(u, v);
    if (vi != none) {
        outEdgeWeightsPerThread[omp_get_thread_num()][u % omp_get_max_threads()][vi] = ew;
        if (!directed && autoCompleteEdges) { // need to adjust both half edges
            const index ui = indexInOutEdgeArrayPerThread(v, u);
            outEdgeWeightsPerThread[omp_get_thread_num()][v % omp_get_max_threads()][ui] = ew;
        }
    } else {
        addHalfEdge(u, v, ew);
    }
}

void GraphBuilder::setInWeight(node u, node v, edgeweight ew) {
    assert(isWeighted());
    assert(isDirected());
    const index vi = indexInInEdgeArrayPerThread(u, v);
    if (vi != none) {
        inEdgeWeightsPerThread[omp_get_thread_num()][u % omp_get_max_threads()][vi] = ew;
    } else {
        addHalfEdge(u, v, ew);
    }
}

void GraphBuilder::increaseOutWeight(node u, node v, edgeweight ew) {
    assert(isWeighted());
    const index vi = indexInOutEdgeArrayPerThread(u, v);
    if (vi != none) {
        outEdgeWeightsPerThread[omp_get_thread_num()][u % omp_get_max_threads()][vi] += ew;
        if (!directed && autoCompleteEdges && u != v) { // need to adjust both half edges
            const index ui = indexInOutEdgeArrayPerThread(v, u);
            outEdgeWeightsPerThread[omp_get_thread_num()][v % omp_get_max_threads()][ui] += ew;
        }
    } else {
        addHalfEdge(u, v, ew);
    }
}

void GraphBuilder::increaseInWeight(node u, node v, edgeweight ew) {
    assert(isWeighted());
    assert(isDirected());
    const index vi = indexInInEdgeArrayPerThread(u, v);
    if (vi != none) {
        inEdgeWeightsPerThread[omp_get_thread_num()][u % omp_get_max_threads()][vi] += ew;
    } else {
        addHalfEdge(u, v, ew);
    }
}

Graph GraphBuilder::completeGraph() {
    Graph G(n, weighted, directed);
#ifdef NETWORKIT_SANITY_CHECKS
    assert(G.checkConsistency());
#endif
    // copy edges and weights
    addHalfEdgesToGraph(G);

    G.setEdgeCount(unsafe, numberOfEdges(G));
    G.setNumberOfSelfLoops(unsafe, selfloops);
    assert(n == G.upperNodeIdBound());
#ifdef NETWORKIT_SANITY_CHECKS
    assert(G.checkConsistency());
#endif
    G.shrinkToFit();

    reset();

    return G;
}

/*void GraphBuilder::toGraphParallel(Graph &G) {
    // basic idea of the parallelization:
    // 1) each threads collects its own data
    // 2) each node collects all its data from all threads

    int maxThreads = omp_get_max_threads();

    using Adjacencylists = std::vector<std::vector<node>>;
    using Weightlists = std::vector<std::vector<edgeweight>>;

    std::vector<Adjacencylists> inEdgesPerThread(maxThreads, Adjacencylists(n));
    std::vector<Weightlists> inWeightsPerThread(weighted ? maxThreads : 0, Weightlists(n));
    std::vector<count> numberOfSelfLoopsPerThread(maxThreads, 0);

    // step 1
    parallelForNodes([&](node v) {
        int tid = omp_get_thread_num();
        for (index i = 0; i < outEdges[v].size(); i++) {
            node u = outEdges[v][i];
            // self loops don't need to be added twice in undirected graphs
            if (directed || u != v) {
                inEdgesPerThread[tid][u].push_back(v);
                if (weighted) {
                    G.addPartialEdge(unsafe, v, u, outEdgeWeights[v][i]);
                    inWeightsPerThread[tid][u].push_back(outEdgeWeights[v][i]);
                } else {
                    G.addPartialEdge(unsafe, v, u);
                }
            }
            if (u == v) {
                numberOfSelfLoopsPerThread[tid]++;
                if (!directed) {
                    inEdgesPerThread[tid][u].push_back(v);
                    if (weighted)
                        inWeightsPerThread[tid][u].push_back(outEdgeWeights[v][i]);
                }
            }
        }
    });

    // step 2
    parallelForNodes([&](node v) {
        count inDeg = 0;
        count outDeg = G.degreeOut(v);
        for (int tid = 0; tid < maxThreads; tid++) {
            inDeg += inEdgesPerThread[tid][v].size();
        }

        assert(inDeg <= n);
        assert(outDeg <= n);
        // allocate enough memory for all edges/weights
        if (directed) {
            G.preallocateDirectedInEdges(v, inDeg);
        } else {
            G.preallocateUndirected(v, outDeg + inDeg);
        }

        // collect 'second' half of the edges
        if (directed) {
            if (weighted) {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++) {
                        G.addPartialInEdge(unsafe, v, inEdgesPerThread[tid][v][j],
                                           inWeightsPerThread[tid][v][j]);
                    }
                }
            } else {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++) {
                        G.addPartialInEdge(unsafe, v, inEdgesPerThread[tid][v][j]);
                    }
                }
            }
        } else {
            if (weighted) {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++) {
                        G.addPartialEdge(unsafe, v, inEdgesPerThread[tid][v][j],
                                         inWeightsPerThread[tid][v][j]);
                    }
                }
            } else {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++) {
                        G.addPartialEdge(unsafe, v, inEdgesPerThread[tid][v][j]);
                    }
                }
            }
        }
    });

    count numSelfLoops = 0;
#pragma omp parallel for reduction(+ : numSelfLoops)
    for (omp_index i = 0; i < static_cast<omp_index>(maxThreads); ++i)
        numSelfLoops += numberOfSelfLoopsPerThread[i];
    G.setNumberOfSelfLoops(unsafe, numSelfLoops);
}

void GraphBuilder::toGraphSequential(Graph &G) {
    G.removeAllEdges();

    std::vector<count> missingEdgesCounts(n, 0);
    count numberOfSelfLoops = 0;

    // 'first' half of the edges
    for (index u = 0; u < outEdges.size(); u++) {
        for (index j = 0; j < outEdges[u].size(); j++) {
            if (weighted) {
                G.addPartialEdge(unsafe, u, outEdges[u][j], outEdgeWeights[u][j]);
            } else {
                G.addPartialEdge(unsafe, u, outEdges[u][j]);
            }
        }
    }

    std::vector<count> outDeg(G.upperNodeIdBound());
    parallelForNodes([&](node v) { outDeg[v] = G.degreeOut(v); });

    // count missing edges for each node
    G.forNodes([&](node v) {
        // increase count of incoming edges for all neighbors
        G.forNeighborsOf(v, [&](node u) {
            if (directed || u != v) {
                ++missingEdgesCounts[u];
            }
            if (u == v) {
                // self loops don't need to be added again
                // but we need to count them
                ++numberOfSelfLoops;
            }
        });
    });

    // 'second' half the edges
    if (directed) {
        // directed: outEdges is complete, missing half edges are the inEdges
        // missingEdgesCounts are our inDegrees
        std::vector<count> inDeg = missingEdgesCounts;

        // reserve the exact amount of space needed
        G.forNodes([&](node v) { G.preallocateDirectedInEdges(v, inDeg[v]); });

        // copy values
        G.forNodes([&](node v) {
            for (index i = 0; i < outDeg[v]; i++) {
                auto nodeWithWeight = G.getIthNeighborWithWeight(unsafe, v, i);
                node u = nodeWithWeight.first;
                if (weighted) {
                    edgeweight ew = nodeWithWeight.second;
                    G.addPartialInEdge(unsafe, u, v, ew);
                } else {
                    G.addPartialInEdge(unsafe, u, v);
                }
            }
        });
    } else {
        // undirected: so far each edge is just saved at one node
        // add it to the other node as well

        // reserve the exact amount of space needed
        G.forNodes([&](node v) { G.preallocateUndirected(v, outDeg[v] + missingEdgesCounts[v]); });

        // copy values
        G.forNodes([&](node v) {
            // the first outDeg[v] edges in G.outEdges[v] are the first half edges
            // we are adding after outDeg[v]
            for (index i = 0; i < outDeg[v]; i++) {
                auto nodeWithWeight = G.getIthNeighborWithWeight(unsafe, v, i);
                node u = nodeWithWeight.first;
                if (u != v) {
                    if (weighted) {
                        edgeweight ew = nodeWithWeight.second;
                        G.addPartialEdge(unsafe, u, v, ew);
                    } else {
                        G.addPartialEdge(unsafe, u, v);
                    }
                } else {
                    // ignore self loops here
                }
            }
        });
    }

    G.setNumberOfSelfLoops(unsafe, numberOfSelfLoops);
}*/

count GraphBuilder::numberOfEdges(const Graph &G) {
    count m = 0;
#pragma omp parallel for reduction(+ : m) if (n > (1 << 20))
    for (omp_index v = 0; v < static_cast<omp_index>(G.upperNodeIdBound()); v++) {
        m += G.degree(v);
    }
    if (G.isDirected()) {
        return m;
    } else {
        // self loops are just counted once
        return (m - selfloops) / 2 + selfloops;
    }
}

void GraphBuilder::addHalfEdgesToGraph(Graph &G) {
    int max_threads = omp_get_max_threads();
#pragma omp parallel num_threads(max_threads)
    {
        int thread_num = omp_get_thread_num();
        std::vector<count> edgeCounts(n / max_threads + 1);
        for (auto &edgesfromThread : outEdgesPerThread) {
            auto &edges = edgesfromThread[thread_num];
            for (HalfEdge edge : edges) {
                ++edgeCounts[edge.source / max_threads];
            }
        }
        for (count i = 0;; ++i) {
            node v = i * max_threads + thread_num;
            if (v >= n)
                break;
            if (!directed) {
                if (!autoCompleteEdges) {
                    // ToDo: allocate tight
                    G.preallocateUndirected(v, (edgeCounts[i] + G.degreeOut(v)) * 2);
                } else {
                    G.preallocateUndirected(v, edgeCounts[i] + G.degreeOut(v));
                }
            } else {
                G.preallocateDirected(v, edgeCounts[i] + G.degreeOut(v),
                                      edgeCounts[i] + G.degreeIn(v));
            }
        }
        for (index i = 0; i < outEdgesPerThread.size(); ++i) {
            auto &edges = outEdgesPerThread[i][thread_num];
            if (weighted) {
                auto &weights = outEdgeWeightsPerThread[i][thread_num];
                for (index j = 0; j < edges.size(); ++j) {
                    G.addPartialOutEdge(Unsafe{}, edges[j].source, edges[j].destination,
                                        weights[j]);
                    if (!autoCompleteEdges) {
                        inEdgesPerThread[edges[j].destination % max_threads][thread_num]
                            .emplace_back(edges[j].destination, edges[j].source);
                        inEdgeWeightsPerThread[edges[j].destination % max_threads][thread_num]
                            .emplace_back(weights[j]);
                    }
                }
            } else {
                for (HalfEdge edge : edges) {
                    G.addPartialOutEdge(Unsafe{}, edge.source, edge.destination);
                    if (!autoCompleteEdges) {
                        inEdgesPerThread[edge.destination % max_threads][thread_num].emplace_back(
                            edge.destination, edge.source);
                    }
                }
            }
        }
        if (!autoCompleteEdges) {
#pragma omp barrier // this is required as inEdgesPerThreads are potentially being added
        }
        if (directed || !autoCompleteEdges) {
            if (directed) {
                std::vector<count> edgeCounts(n / max_threads + 1);
                auto &edgesfromThread = inEdgesPerThread[thread_num];
                for (auto &edges : edgesfromThread) {
                    for (HalfEdge edge : edges) {
                        ++edgeCounts[edge.source / max_threads];
                    }
                }
                for (count i = 0;; ++i) {
                    node v = i * max_threads + thread_num;
                    if (v >= n)
                        break;
                    G.preallocateDirectedInEdges(v, edgeCounts[i] + G.degreeIn(v));
                }
                for (auto &edges : edgesfromThread) {
                    for (HalfEdge edge : edges) {
                        G.addPartialInEdge(Unsafe{}, edge.source, edge.destination);
                    }
                }
            } else { // collect "second" half of the edges
                auto &edgesfromThread = inEdgesPerThread[thread_num];
                for (index i = 0; i < edgesfromThread.size(); ++i) {
                    auto &edges = edgesfromThread[i];
                    if (weighted) {
                        auto &weights = inEdgeWeightsPerThread[thread_num][i];
                        for (index j = 0; j < edges.size(); ++j) {
                            if (edges[j].source != edges[j].destination)
                                G.addPartialOutEdge(Unsafe{}, edges[j].source, edges[j].destination,
                                                    weights[j]);
                        }
                    } else {
                        for (HalfEdge edge : edges) {
                            if (edge.source != edge.destination)
                                G.addPartialOutEdge(Unsafe{}, edge.source, edge.destination);
                        }
                    }
                }
            }
        }
        if (!autoCompleteEdges) {
            G.setNumberOfSelfLoops(unsafe, selfloops);
        } else {
            G.setNumberOfSelfLoops(unsafe, selfloops / 2);
        }
    }
}

} /* namespace NetworKit */
