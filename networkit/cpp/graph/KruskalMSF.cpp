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

struct MyEdge : WeightedEdge {
    using WeightedEdge::WeightedEdge; // Inherit constructors from WeightedEdge

    bool operator<(const MyEdge &other) const noexcept {
        return this->weight > other.weight; // Note the switch in the operator!
    }
};

KruskalMSF::KruskalMSF(const Graph &G) : SpanningForest(G) {}

void KruskalMSF::run() {
    if (G->isWeighted()) {
        count z = G->upperNodeIdBound();
        forest = GraphTools::copyNodes(*G);
        UnionFind uf(z);

        // sort edges in decreasing weight order
        std::vector<MyEdge> sortedEdges(G->numberOfEdges());
        std::transform(G->edgeWeightRange().begin(), G->edgeWeightRange().end(),
                       sortedEdges.begin(), [](const auto &edge) -> MyEdge {
                           return {edge.u, edge.v, edge.weight};
                       });
        Aux::Parallel::sort(sortedEdges.begin(), sortedEdges.end());

        // process in decreasing weight order
        for (const auto &e : sortedEdges) {
            DEBUG("process edge (", e.u, ", ", e.v, ") with weight ", e.weight);

            // if edge does not close cycle, add it to tree
            if (uf.find(e.u) != uf.find(e.v)) {
                forest.addEdge(e.u, e.v);
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
