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
class DynBetweenness: public Centrality {

public:
    /**
     * Constructs the Betweenness class for the given Graph @a G.
     *
     * @param G The graph.
     */
    DynBetweenness(const Graph& G, bool storePredecessors = false);

    /**
     * Runs the static betweenness centrality algorithm on the initial graph.
     */
    void run() override;

    /**
    * Updates the betweenness centralities after a batch of edge insertions on the graph.
    *
    * @param batch The batch of edge insertions.
    */
    void update(GraphEvent e);


    edgeweight distance(node s, node t) const;

    double dependency(node s, node t) const;

    count nPaths(node s, node t) const;



protected:
    void updateWeighted(GraphEvent e);
    void updateUnweighted(GraphEvent e);
    std::vector<count> maxDistance;
    std::vector<std::vector<count>> npaths;
    std::vector<std::vector<edgeweight>> distances;
    std::vector<std::vector<double>> dependencies;
    bool storePreds = false;
    std::vector<std::vector<std::vector<node>>> predecessors;
};

inline edgeweight DynBetweenness::distance(node s, node t) const {
    return distances[s][t];
}

inline double DynBetweenness::dependency(node s, node t) const {
    return dependencies[s][t];
}

inline count DynBetweenness::nPaths(node s, node t) const {
    return npaths[s][t];
}

} /* namespace NetworKit */

#endif /* DYNBETW_H_ */
