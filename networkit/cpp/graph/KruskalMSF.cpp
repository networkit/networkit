/*
 * KruskalMSF.cpp
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/KruskalMSF.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace NetworKit {

KruskalMSF::KruskalMSF(const Graph &G) : SpanningForest(G) {}

void KruskalMSF::run() {
    if (G->isWeighted()) {
        count z = G->upperNodeIdBound();
        forest = GraphTools::copyNodes(*G);
        UnionFind uf(z);

        // sort edges in non-decreasing weight order
        std::vector<WeightedEdge> sortedEdges(G->numberOfEdges());
        std::copy(G->edgeWeightRange().begin(), G->edgeWeightRange().end(), sortedEdges.begin());
        Aux::Parallel::sort(
            sortedEdges.begin(), sortedEdges.end(),
            [](const WeightedEdge &u, const WeightedEdge &v) { return u.weight < v.weight; });

        // process in non-decreasing weight order
        for (const auto &e : sortedEdges) {
            DEBUG("process edge (", e.u, ", ", e.v, ") with weight ", e.weight);

            // if edge does not close cycle, add it to tree
            if (uf.find(e.u) != uf.find(e.v)) {
                forest.addEdge(e.u, e.v);
                totalWeight += e.weight;
                uf.merge(e.u, e.v);
            }
        }
    } else {
        SpanningForest sf(*G);
        sf.run();
        forest = sf.getForest();
    }

    hasRun = true;
}

} /* namespace NetworKit */
