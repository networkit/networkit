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

void APSP::run() {
    const count n = G.upperNodeIdBound();
    distances.assign(n, std::vector<edgeweight>(n));

    sssps.resize(omp_get_max_threads());
#pragma omp parallel
    {
        omp_index i = omp_get_thread_num();
        if (G.isWeighted())
            sssps[i] = std::unique_ptr<SSSP>(new Dijkstra(G, 0, false));
        else
            sssps[i] = std::unique_ptr<SSSP>(new BFS(G, 0, false));
    }

    G.parallelForNodes([&](node source) {
        auto sssp = sssps[omp_get_thread_num()].get();
        sssp->setSource(source);
        sssp->run();
        distances[source] = sssp->getDistances();
    });

    hasRun = true;
}

} /* namespace NetworKit */
