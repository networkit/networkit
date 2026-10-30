#ifndef NETWORKIT_ISOMORPHISM_PARALLEL_RI_HPP_
#define NETWORKIT_ISOMORPHISM_PARALLEL_RI_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/RI.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {

/**
 * @ingroup isomorphism
 * Parallel version of @ref RI.
 *
 * ParallelRI finds the same matches as @ref RI, but their order may differ from run to run. Under
 * a cap on the number of matches, which matches are kept may differ as well. The number of
 * workers is the global thread count at the start of @ref run(), see `Aux::setNumberOfThreads()`.
 * Busy workers publish their oldest partial mappings in small batches, and idle workers steal them
 * one at a time. A
 * @ref SubgraphIsomorphism::MatchCallback makes the workers wait for each other, and a
 * @ref SubgraphIsomorphism::ParallelMatchCallback avoids this.
 *
 * The implementation is based on
 *
 * Kimmig, R., Meyerhenke, H., & Strash, D. (2017).
 * Shared Memory Parallel Subgraph Enumeration.
 * IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW).
 */
class ParallelRI final : public SubgraphIsomorphism {

public:
    /**
     * @param pattern The graph to look for. Must not contain self-loops.
     * @param target The graph to look in. Must agree with @a pattern on directedness.
     * @param variant Plain RI or RI-DS. See @ref RI::Variant.
     * @param semantics Whether matches must be induced. See @ref SubgraphIsomorphism::Semantics.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     * @param domainFilterThreshold Only for @ref RI::Variant::RI_DS. Decides which domains filter
     * the candidates, see the constructor of @ref RI.
     *
     * @throws std::invalid_argument If @a variant or @a semantics is not a known value, or if
     * @a domainFilterThreshold is not between 0 and 1. @ref run() checks the graphs, since they may
     * change after construction.
     */
    ParallelRI(const Graph &pattern, const Graph &target, RI::Variant variant = RI::Variant::RI,
               Semantics semantics = Semantics::INDUCED, count maxMatches = 0,
               double domainFilterThreshold = 0.8);

    void run() override;

private:
    RI::Variant variant;
    double domainFilterThreshold;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_PARALLEL_RI_HPP_
