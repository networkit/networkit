/*
 * SelectiveCommunityDetector.hpp
 *
 *  Created on: 15.05.2013
 *      Author: cls
 */

#ifndef NETWORKIT_SCD_SELECTIVE_COMMUNITY_DETECTOR_HPP_
#define NETWORKIT_SCD_SELECTIVE_COMMUNITY_DETECTOR_HPP_

#include <map>
#include <set>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class SelectiveCommunityDetector {

public:

    /**
     * @param[in] G the input graph.
     */
    SelectiveCommunityDetector(const Graph& G) : G(&G) {}
    virtual ~SelectiveCommunityDetector() = default;

    /**
     * Detect communities for given seed nodes.
     * @return a mapping from seed node to community (as a set of nodes)
     */
    virtual std::map<node, std::set<node>> run(const std::set<node>& seeds) = 0;

protected:
    const Graph* G;
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_SELECTIVE_COMMUNITY_DETECTOR_HPP_
