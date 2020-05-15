/*
 * APSP.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#include <omp.h>

#include <networkit/distance/APSP.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>

namespace NetworKit {

APSP::APSP(const Graph &G) : Algorithm(), G(G) {}

APSP::APSP(const Graph &G, std::vector<count> &srcs) : Algorithm(), G(G), srcs(srcs) {}

void APSP::run() {
    const count n = G.upperNodeIdBound();
    const bool src_all = srcs.size() == 0;
    const count num_srcs = src_all ? n : srcs.size();
    std::vector<edgeweight> distanceVector(n, 0.0);
    distances.resize(num_srcs, distanceVector);

    count nThreads = omp_get_max_threads();
    sssps.resize(nThreads);
#pragma omp parallel
    {
        omp_index i = omp_get_thread_num();
        if (G.isWeighted())
            sssps[i] = std::unique_ptr<SSSP>(new Dijkstra(G, 0, false));
        else
            sssps[i] = std::unique_ptr<SSSP>(new BFS(G, 0, false));
    }

#pragma omp parallel for schedule(dynamic)
    for (omp_index source_idx = 0; source_idx < num_srcs; ++source_idx) {
        auto sssp = sssps[omp_get_thread_num()].get();
        count source = src_all ? source_idx : srcs[source_idx];
        sssp->setSource(source);
        sssp->run();
        distances[source_idx] = sssp->getDistances();
    }

    hasRun = true;
}

std::string APSP::toString() const {
    return "All-Pairs Shortest Path Algorithm";
}

} /* namespace NetworKit */
