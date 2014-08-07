/*
 * DynSSSP.h
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#ifndef DYNSSSP_H_
#define DYNSSSP_H_

#include <set>

#include "Graph.h"
#include "../dynamics/GraphEvent.h"
#include "SSSP.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Interface for dynamic single-source shortest path algorithms.
 */
class DynSSSP: public SSSP {

public:
    /**
     * The algorithm computes a dynamic SSSP starting from the specified
     * source vertex.
     *
     * @param	graph		input graph.
     * @param	source	    source vertex.
     */
    DynSSSP(const Graph& G, node source);
    /**
    * Updates the betweenness centralities after a batch of edge insertions on the graph.
    *
    * @param batch The batch of edge insertions.
    */
    virtual void update(const std::vector<GraphEvent>& batch) = 0;
    /**
    * Returns true or false depending on whether the node previoulsy specified
    * with setTargetNode has been modified by the udate or not.
    *
    * @param batch The batch of edge insertions.
    */
    bool modified();
    /**
    * Set a target node to be `observed` during the update. If a node t is set as
    * target before the update, the function modified() will return true or false
    * depending on whether node t has been modified by the update.
    *
    * @param t Node to be `observed`.
    */
    void setTargetNode(const node t = 0);

protected:
    bool mod = false;
    node target;
};

inline bool DynSSSP::modified() {
    return mod;
}

inline void DynSSSP::setTargetNode(const node t) {
    target = t;
}

} /* namespace NetworKit */

#endif /* DYNSSSP_H_ */
