/*
 * SpanningForest.hpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#ifndef NETWORKIT_GRAPH_SPANNING_FOREST_HPP_
#define NETWORKIT_GRAPH_SPANNING_FOREST_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Base class for spanning forest/tree algorithms.
 */
template <typename GraphT>
class GenericSpanningForest : public Algorithm {
public:
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

protected:
    const GraphT *G;
    GraphT forest;

    static GraphT copyNodes(const GraphT &G);

public:
    GenericSpanningForest(const GraphT &G) : G(&G) {}

    void run() override;

    /**
     * @return Forest computed by run method.
     * Note: So far no explicit check if run method has been invoked before.
     */
    const GraphT &getForest() {
        assureFinished();
        return forest;
    }
};

using SpanningForest = GenericSpanningForest<Graph>;

} /* namespace NetworKit */

#include <networkit/graph/SpanningForestImpl.hpp>

#endif // NETWORKIT_GRAPH_SPANNING_FOREST_HPP_
