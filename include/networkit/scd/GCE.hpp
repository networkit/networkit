/* GCE.hpp
 *
 * Created on: 06.05.2013
 * Author: cls
 */

#ifndef NETWORKIT_SCD_GCE_HPP_
#define NETWORKIT_SCD_GCE_HPP_

#include <unordered_set>

#include <networkit/auxiliary/SetIntersector.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * The Greedy Community Expansion algorithm.
 *
 * Greedily adds nodes from the shell to improve community quality.
 */
class GCE : public SelectiveCommunityDetector {

public:
    GCE(const Graph &g, std::string objective);

    /**
     * @param[in]  seeds  seed nodes
     *
     * @param[out]        community as a set of nodes
     */
    std::set<node> expandOneCommunity(const std::set<node> &seeds) override;

    // inherit method from parent class.
    using SelectiveCommunityDetector::expandOneCommunity;

private:
    std::string objective; // name of objective function
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_GCE_HPP_
