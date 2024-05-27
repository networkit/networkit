/*
 * Hypergraph.cpp
 *
 *  Created on: 24.05.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

Hypergraph::Hypergraph(count n, count m, bool weighted, bool directed)
    : n(n), m(m), z(n), omega(0),

      weighted(weighted), // indicates whether the graph is weighted or not
      directed(directed), // indicates whether the graph is directed or not

      nodeExists(n, true), edgeExists(m, true),

      nodeWeights(weighted ? n : 0), edgeWeights(weighted ? m : 0),

      nodeAttributeMap(this), edgeAttributeMap(this),

      /* for directed graphs nodeHeadIncidence stores an incidence list only considering
       * edges where the node is present as head, for undirected hypergraphs nodeHeadIncidence is
       * not used. The same is true for edgeHeadIncidence
       */
      nodeHeadIncidence(directed ? n : 0), edgeHeadIncidence(directed ? m : 0) {}
} // namespace NetworKit
