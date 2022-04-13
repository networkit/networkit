/*
 * SPSP.cpp
 *
 *  Created on: 23.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <memory>
#include <omp.h>

#include <networkit/distance/BFS.hpp>
#include <networkit/distance/BidirectionalBFS.hpp>
#include <networkit/distance/BidirectionalDijkstra.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/MultiTargetBFS.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/distance/SPSP.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

void SPSP::run() {
    distances.resize(sources.size());

    if (targets.empty())
        runWithoutTargets();
    else
        runWithTargets();

    hasRun = true;
}

void SPSP::runWithoutTargets() {
#pragma omp parallel
    {
        std::unique_ptr<SSSP> sssp;
        if (G->isWeighted())
            sssp = std::unique_ptr<SSSP>(new Dijkstra(*G, 0, false));
        else
            sssp = std::unique_ptr<SSSP>(new BFS(*G, 0, false));

#pragma omp for
        for (omp_index i = 0; i < sources.size(); ++i) {
            sssp->setSource(sources[i]);
            sssp->run();
            distances[i] = sssp->getDistances();
        }
    }
}

void SPSP::runWithTargets() {
#pragma omp parallel
    {
        std::unique_ptr<STSP> stsp;

        {
            // Init algos with a dummy source node, we have to change source anyway later
            node dummy = *G->nodeRange().begin();
            if (G->isWeighted())
                stsp = std::unique_ptr<STSP>(
                    new MultiTargetDijkstra(*G, dummy, targets.begin(), targets.end()));
            else
                stsp = std::unique_ptr<STSP>(
                    new MultiTargetBFS(*G, dummy, targets.begin(), targets.end()));
        }

#pragma omp for
        for (omp_index i = 0; i < static_cast<omp_index>(sources.size()); ++i) {

            stsp->setSource(sources[i]);
            stsp->run();

            const auto &dists = stsp->getDistances();
            const auto &stspTgtIdx = stsp->getTargetIndexMap();
            distances[i].reserve(targets.size());

            for (node target : targets)
                distances[i].push_back(dists[stspTgtIdx.at(target)]);
        }
    }
}

} // namespace NetworKit
