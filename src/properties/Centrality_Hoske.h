/*
 * Centrality_Hoske.h
 *
 *  Created on: 05.12.2013
 *      Author: dhoske
 */

#ifndef CENTRALITYHOSKE_H_
#define CENTRALITYHOSKE_H_

#include "../graph/Graph.h"

namespace NetworKit {

/** Computes the betweenness centrality of the nodes in G. */
std::vector<double> betweennessCentrality_Hoske(const Graph& G);

} /* namespace NetworKit */
#endif /* CENTRALITYHOSKE_H_ */
