
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

    void run() override;

    /**
     * @return Forest computed by run method.
     * Note: So far no explicit check if run method has been invoked before.
     */
    const Graph &getForest() {
        assureFinished();
        return forest;
    }

};
}
#endif //NETWORKIT_GRAPH_ORIM_MST_H
