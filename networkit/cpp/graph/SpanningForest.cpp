/*
 * SpanningForest.cpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#include "SpanningForest.h"

namespace NetworKit {

SpanningForest::SpanningForest(const Graph& G): G(G) {

}

Graph SpanningForest::getForest() {
	return forest;
}

} /* namespace NetworKit */
