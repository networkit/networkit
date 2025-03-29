
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

class PrimMST : public SpanningForest {

    PrimMST(const Graph &G) : SpanningForest(G) {}

    void run() override;

    const Graph &getForest() {
        assureFinished();
        return forest;
    }
private:
    static constexpr edgeweight infiniteWeight = std::numeric_limits<edgeweight>::max();
};
}
#endif //NETWORKIT_GRAPH_ORIM_MST_H
