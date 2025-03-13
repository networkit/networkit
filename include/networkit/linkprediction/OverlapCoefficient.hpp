/*
 * OverlapCoefficient.hpp
 *
 *  Created on: 13.03.2025
 *      Authors: Fabian Brandt-Tumescheit
 */

#ifndef NETWORKIT_LINKPREDICTION_OVERLAP_COEFFICIENT_HPP_
#define NETWORKIT_LINKPREDICTION_OVERLAP_COEFFICIENT_HPP_

#include <cmath>

#include <networkit/linkprediction/LinkPredictor.hpp>
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Assigns a distance value to pairs of nodes according to the
 * overlap coefficient of their neighborhoods.
 */
class OverlapCoefficient final : public LinkPredictor {
    /**
     * Returns the Neighborhood Distance index for the given node-pair (@a u, @a v).
     * @param u First node
     * @param v Second node
     * @return the Neighborhood Distance index for the given node-pair (@a u, @a v)
     */
    double runImpl(node u, node v) override {
        count uNeighborhood = G->degree(u);
        count vNeighborhood = G->degree(v);
        count intersection = NeighborhoodUtility::getCommonNeighbors(*G, u, v).size();
        return ((double)intersection) / (std::min(uNeighborhood, vNeighborhood));
    }

public:
    using LinkPredictor::LinkPredictor;
};

} /* namespace NetworKit */
#endif // NETWORKIT_LINKPREDICTION_OVERLAP_COEFFICIENT_HPP_
