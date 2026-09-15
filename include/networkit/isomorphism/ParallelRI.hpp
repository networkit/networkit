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
 * See @ref SubgraphIsomorphism for the definition of a match and how to use the class, and
 * @ref RI for the algorithm. ParallelRI finds the same matches as @ref RI, but the order of the
 * matches may differ from run to run.
 *
 * Every worker expands partial mappings from its own queue, depth first. To balance the load, a
 * worker publishes batches of its oldest partial mappings, and an idle worker steals from the
 * batches of a randomly chosen other worker. A token that is passed around a ring of workers
 * detects when all of them are idle. The number of workers is the global thread count, see
 * `Aux::setNumberOfThreads()`.
 *
 * A @ref SubgraphIsomorphism::MatchCallback is never called concurrently, so the workers wait for
 * each other to call it. A @ref SubgraphIsomorphism::ParallelMatchCallback avoids this.
 *
 * The search can be interrupted with CTRL+C.
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
     */
    ParallelRI(const Graph &pattern, const Graph &target, RI::Variant variant = RI::Variant::RI,
               Semantics semantics = Semantics::INDUCED, count maxMatches = 0);

    /**
     * Runs the search. Query the result with @ref getMatches(), @ref numberOfMatches() or
     * @ref hasMatch().
     */
    void run() override;

    /**
     * Returns the number of workers @ref run() uses, which is the global thread count.
     *
     * @return the number of workers.
     */
    count numberOfWorkers() const override;

private:
    RI::Variant variant;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_PARALLEL_RI_HPP_
