/*
 * DynApproxBetweenness.h
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

#ifndef DYNAPPROXBETW_H_
#define DYNAPPROXBETW_H_

#include "Centrality.h"
#include "DynCentrality.h"
#include "../dynamics/GraphEvent.h"
#include "../graph/DynSSSP.h"

#include <math.h>
#include <algorithm>
#include <memory>
#include <omp.h>

namespace NetworKit {

/**
 * @ingroup graph
 * Interface for dynamic approximated betweenness centrality algorithm.
 */
class DynApproxBetweenness: public Centrality, public DynCentrality {

public:
    /**
      * The algorithm approximates the betweenness of all vertices so that the scores are
      * within an additive error @a epsilon with probability at least (1- @a delta).
      * The values are normalized by default.
      *
      * @param	G			the graph
      * @param  storePredecessors   keep track of the lists of predecessors?
      * @param	epsilon		maximum additive error
      * @param	delta		probability that the values are within the error guarantee
     */
    DynApproxBetweenness(const Graph& G, double epsilon=0.01, double delta=0.1, bool storePredecessors = true);

    /**
     * Runs the static approximated betweenness centrality algorithm on the initial graph.
     */
    void run() override;

    /**
    * Updates the betweenness centralities after a batch of edge insertions on the graph.
    *
    * @param batch The batch of edge insertions.
    */
    void update(const std::vector<GraphEvent>& batch);

    /**
    * Get number of path samples used for last calculation
    */
    count getNumberOfSamples();

private:
    bool storePreds = true;
    double epsilon; //!< maximum error
    double delta;
    count r;
    std::vector<std::unique_ptr<DynSSSP>> sssp;
    std::vector<node> u;
    std::vector<node> v;
    std::vector <std::vector<node>> sampledPaths;
};

} /* namespace NetworKit */

#endif /* DYNAPPROXBETW_H_ */
