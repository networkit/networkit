/*
 * DynBetweenness.h
 *
 *  Created on: 29.07.2014
 *      Author: ebergamini
 */

#ifndef DYNBETW_H_
#define DYNBETW_H_

#include "Centrality.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Interface for dynamic betweenness centrality algorithms.
 */
class DynBetweenness: public NetworKit::Centrality {

public:
    /**
     * Constructs the Betweenness class for the given Graph @a G. If the betweenness scores should be normalized,
     * then set @a normalized to <code>true</code>.
     *
     * @param G The graph.
     * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
     */
    DynBetweenness(const Graph& G);

    /**
     * Runs the static betweenness centrality algorithm on the initial graph.
     */
    void run() override;

    /**
    * Updates the betweenness centralities after a batch of edge insertions on the graph.
    *
    * @param batch The batch of edge insertions.
    */
    void update(const GraphEvent e);

protected:
        std::vector<count> maxDistance;
        std::vector<std::vector<count>> npaths;
        std::vector<std::vector<edgeweight>> distances;
        std::vector<std::vector<double>> dependencies;
};

} /* namespace NetworKit */

#endif /* DYNBETW_H_ */
