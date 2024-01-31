
#include <omp.h>
#include <stdexcept>

#include <networkit/graph/ParallelGraphBuilder.hpp>

namespace NetworKit {

ParallelGraphBuilder::ParallelGraphBuilder(count n, bool directed)
    : n(n), directed(directed), outEdgesPerThread(omp_get_max_threads()),
      inEdgesPerThread(directed ? omp_get_max_threads() : 0) {
    for (int i = 0; i < omp_get_max_threads(); i++) {
        outEdgesPerThread[i].resize(omp_get_max_threads());
    }
    if (directed) {
        for (int i = 0; i < omp_get_max_threads(); i++) {
            inEdgesPerThread[i].resize(omp_get_max_threads());
        }
    }
}

void ParallelGraphBuilder::addEdgeParallel(index a, index b) {
    if (directed) {
        addHalfOutEdgeParallel(Unsafe{}, a, b);
        addHalfInEdgeParallel(Unsafe{}, b, a);
    } else {
        addHalfEdgeParallel(Unsafe{}, a, b);
        addHalfEdgeParallel(Unsafe{}, b, a);
    }
}

void ParallelGraphBuilder::addHalfEdgeParallel(Unsafe, index a, index b) {
    assert(a != b);
    outEdgesPerThread[omp_get_thread_num()][a % omp_get_max_threads()].emplace_back(HalfEdge{a, b});
}

void ParallelGraphBuilder::addHalfInEdgeParallel(Unsafe, index a, index b) {
    assert(a != b);
    inEdgesPerThread[omp_get_thread_num()][a % omp_get_max_threads()].emplace_back(HalfEdge{a, b});
}

void ParallelGraphBuilder::addHalfOutEdgeParallel(Unsafe, index a, index b) {
    assert(a != b);
    outEdgesPerThread[omp_get_thread_num()][a % omp_get_max_threads()].emplace_back(HalfEdge{a, b});
}

void ParallelGraphBuilder::addEdgesToGraphParallel(Unsafe, Graph &G) {
    count numberOfHalfEdges = 2 * G.numberOfEdges();
#pragma omp parallel num_threads(omp_get_max_threads())
    {
        std::vector<count> edgeCounts(n / omp_get_max_threads() + 1);
        for (auto &edgesfromThread : outEdgesPerThread) {
            auto &edges = edgesfromThread[omp_get_thread_num()];
            for (auto edge : edges) {
                edgeCounts[edge.source / omp_get_max_threads()]++;
            }
        }
        for (count i = 0;; i++) {
            node v = i * omp_get_max_threads() + omp_get_thread_num();
            if (v >= n)
                break;
            G.preallocateUndirected(v, edgeCounts[i] + G.degreeOut(v));
        }
        for (auto &edgesfromThread : outEdgesPerThread) {
            auto &edges = edgesfromThread[omp_get_thread_num()];
            for (auto edge : edges) {
                G.addPartialOutEdge(Unsafe{}, edge.source, edge.destination);
            }
#pragma omp atomic
            numberOfHalfEdges += edges.size();
        }
        if (directed) {
            for (count &edgeCount : edgeCounts) {
                edgeCount = 0;
            }
            for (auto &edgesfromThread : inEdgesPerThread) {
                auto &edges = edgesfromThread[omp_get_thread_num()];
                for (auto edge : edges) {
                    edgeCounts[edge.source / omp_get_max_threads()]++;
                }
            }
            for (count i = 0;; i++) {
                node v = i * omp_get_max_threads() + omp_get_thread_num();
                if (v >= n)
                    break;
                G.preallocateDirectedInEdges(v, edgeCounts[i] + G.degreeIn(v));
            }
            for (auto &edgesfromThread : inEdgesPerThread) {
                auto &edges = edgesfromThread[omp_get_thread_num()];
                for (auto edge : edges) {
                    G.addPartialInEdge(Unsafe{}, edge.source, edge.destination);
                }
#pragma omp atomic
                numberOfHalfEdges += edges.size();
            }
        }
    }
    G.setEdgeCount(unsafe, numberOfHalfEdges / 2);
}

Graph ParallelGraphBuilder::toGraphParallel() {
    Graph G(n, false, directed);
    addEdgesToGraphParallel(Unsafe{}, G);
    return G;
}

} /* namespace NetworKit */