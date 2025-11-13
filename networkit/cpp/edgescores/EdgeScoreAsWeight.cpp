/*
 * EdgeScoreAsWeight.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include "networkit/edgescores/EdgeScore.hpp"
#include <networkit/edgescores/EdgeScoreAsWeight.hpp>

namespace NetworKit {

EdgeScoreAsWeight::EdgeScoreAsWeight(const Graph &G, const std::vector<double> &score, bool squared,
                                     edgeweight offset, edgeweight factor)
    : G(&G), score(&score), squared(squared), offset(offset), factor(factor) {}

Graph EdgeScoreAsWeight::calculate() {
    WARN("The class EdgeScoreAsWeight is deprecated and will be removed in future releases. Use "
         "EdgeScore<T>::calculate(...) instead.");
    return EdgeScore<edgeweight>(*G).calculate(squared, offset, factor);
}

} // namespace NetworKit
