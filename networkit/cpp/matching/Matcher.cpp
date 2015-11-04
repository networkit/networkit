/*
 * Matcher.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Matcher.h"

namespace NetworKit {

Matcher::Matcher(const Graph& G): G(G), M(G.upperNodeIdBound()) {

}

Matching Matcher::getMatching() const {
    return M;
}



} /* namespace NetworKit */
