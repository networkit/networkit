
/*  PrimMST.hpp
 *
 *	Created on: 29.03.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#ifndef NETWORKIT_GRAPH_ORIM_MST_H
#define NETWORKIT_GRAPH_ORIM_MST_H
#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/SpanningForest.hpp>

namespace NetworKit {

class PrimMSF : public SpanningForest {
public:
    PrimMSF(const Graph &G) : SpanningForest(G) {}

    void run() override;

    edgeweight getTotalWeight() const {
        assureFinished();
        if (G->isWeighted())
            return totalWeight;
        return static_cast<edgeweight>(forest.numberOfEdges());
    }
private:
    static constexpr edgeweight infiniteWeight = std::numeric_limits<edgeweight>::max();
    edgeweight totalWeight = 0;
};
}
#endif //NETWORKIT_GRAPH_ORIM_MST_H
