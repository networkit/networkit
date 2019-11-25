/*
 * ParallelMatcher.h
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_
#define NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_

#include <set>
#include <algorithm>

#include <networkit/matching/Matcher.hpp>

namespace NetworKit {


/**
 * @ingroup matching
 * LocalMax matching similar to the one described in the EuroPar13 paper
 * by the Sanders group (Birn, Osipov, Sanders, Schulz, Sitchinava)
 */
class LocalMaxMatcher: public Matcher {
public:

    LocalMaxMatcher(const Graph& G);


    virtual void run();

protected:

};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_
