/*
 * QualityMeasure.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_COMMUNITY_QUALITY_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_QUALITY_MEASURE_HPP_

#include <networkit/structures/Partition.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for all clustering quality measures.
 */
class QualityMeasure {


public:
    virtual double getQuality(const Partition& zeta, const Graph& G) = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_QUALITY_MEASURE_HPP_
