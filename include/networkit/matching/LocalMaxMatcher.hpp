/*
 * LocalMaxMatcher.hpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_
#define NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_

#include <algorithm>
#include <set>

#include <networkit/matching/Matcher.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * LocalMax matching similar to the one described in the EuroPar13 paper
 * by the Sanders group (Birn, Osipov, Sanders, Schulz, Sitchinava)
 */
class LocalMaxMatcher final : public Matcher {
public:

    LocalMaxMatcher(const Graph& G);

    void run() override;

};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_
