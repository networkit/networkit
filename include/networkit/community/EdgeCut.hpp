/*
 * EdgeCut.h
 *
 *  Created on: Jun 20, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_COMMUNITY_EDGE_CUT_HPP_
#define NETWORKIT_COMMUNITY_EDGE_CUT_HPP_

#include <networkit/community/QualityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 */
class EdgeCut: public QualityMeasure {
public:
    virtual double getQuality(const Partition& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_EDGE_CUT_HPP_
