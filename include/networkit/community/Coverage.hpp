/*
 * Coverage.h
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_COMMUNITY_COVERAGE_HPP_
#define NETWORKIT_COMMUNITY_COVERAGE_HPP_

#include <networkit/community/QualityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Coverage is the fraction of intra-cluster edges.
 */
class Coverage: public QualityMeasure {
public:

    virtual double getQuality(const Partition& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_COVERAGE_HPP_
