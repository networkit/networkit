#ifndef NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_

// Private header of the isomorphism module. Not installed, not part of the public API. Shared by
// RI.cpp and ParallelRI.cpp.

#include <string>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/isomorphism/RI.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

#include "MatchReporter.hpp"
#include "SearchGraph.hpp"

namespace NetworKit {
namespace IsomorphismDetails {

/**
 * The RI search, shared by @ref RI and @ref ParallelRI.
 *
 * The preprocessing, @ref computeDomains() followed by @ref computeOrdering(), depends only on the
 * inputs. The driver computes it once, and all workers of @ref ParallelRI share the read-only
 * result. The search runs in one of two ways. @ref run() enumerates all matches recursively and is
 * used by @ref RI. @ref expand() produces the children of one @ref State and is used by the
 * workers of @ref ParallelRI, so that unexplored states can be stolen. Both use
 * @ref candidatesFor(), @ref consistent() and @ref reportMapping(), so they find the same matches.
 */
class RIImpl {

public:
    /**
     * The fixed order in which pattern nodes are mapped. `order[i]` is the pattern node mapped at
     * depth @a i. `parent[i]` is the position in `order` of an earlier neighbour of `order[i]`, or
     * @ref none if position @a i starts a new connected component. The candidates for position
     * @a i are drawn from the neighbourhood of the parent's image.
     */
    struct Ordering {
        std::vector<node> order;
        std::vector<index> parent;
    };

    /**
     * RI-DS candidate domains, indexed by pattern node id. Empty under plain RI.
     */
    struct Domains {
        /// `ofPatternNode[pu]` holds the target nodes that @a pu can be mapped to, in ascending
        /// order.
        std::vector<std::vector<node>> ofPatternNode;

        /// Whether the refinement removed enough of a domain to make intersecting target
        /// neighbourhoods with it worthwhile. See MinSweepYieldForSliceIntersection in RIImpl.cpp.
        std::vector<bool> earnsItsKeep;

        /// True if some domain is empty, so that there are no matches at all.
        bool anyEmpty = false;
    };

    /**
     * Computes the RI-DS domains. Returns empty domains under plain RI.
     *
     * A target node enters the domain of pattern node @a pu if it exists, its degrees are at least
     * those of @a pu, and its label is compatible. One refinement pass then keeps a candidate only
     * if, for every pattern arc at @a pu, it has a target arc of the same direction and a
     * compatible label into the domain of the other endpoint. Finally, forward checking removes
     * the target node of every single-element domain from all other domains, and repeats this for
     * domains that become single-element.
     *
     * @param pattern Snapshot of the pattern.
     * @param target Snapshot of the target.
     * @param patternNodeLabels Empty if the search is unlabelled.
     * @param targetNodeLabels Empty if the search is unlabelled.
     * @param variant The RI variant to run.
     * @return the domains.
     */
    static Domains computeDomains(const SearchGraph &pattern, const SearchGraph &target,
                                  const std::vector<index> &patternNodeLabels,
                                  const std::vector<index> &targetNodeLabels, RI::Variant variant);

    /**
     * A partial mapping, that is, a node of the search tree. @ref ParallelRI steals these between
     * workers.
     *
     * `mapping[i]` is the target node that `order[i]` is mapped to. It has one entry per position
     * in the order, with @ref none beyond `depth`. `depth` is the number of mapped positions, and
     * `nextCandidate` is the index in the candidate list at which @ref expand() resumes.
     */
    struct State {
        std::vector<node> mapping;
        count depth = 0;
        index nextCandidate = 0;
    };

    /**
     * Computes the order in which pattern nodes are mapped.
     *
     * The order starts at a node of maximum degree and repeatedly appends the unordered node that
     * maximizes three counts in lexicographic order: the arcs into the ordered nodes, the ordered
     * nodes reachable in two hops through an unordered node, and the unordered neighbours that are
     * adjacent to no ordered node. Remaining ties go to the smallest node id.
     *
     * Under RI-DS, pattern nodes with a single-element domain are ordered first, and a tie on all
     * three counts goes to the smaller domain before the node id.
     *
     * @param pattern Snapshot of the pattern.
     * @param domains Result of @ref computeDomains(); empty under plain RI.
     * @return the matching order.
     */
    static Ordering computeOrdering(const SearchGraph &pattern, const Domains &domains);

    /**
     * @return true if the pattern has more nodes or a larger maximum degree than the target, in
     * which case there are no matches.
     */
    static bool patternCannotFit(const SearchGraph &pattern, const SearchGraph &target);

    /**
     * @param pattern Snapshot of the pattern, built with the adjacency matrix.
     * @param target Snapshot of the target, built without it.
     * @param patternNodeLabels Empty if the search is unlabelled.
     * @param targetNodeLabels Empty if the search is unlabelled.
     * @param ordering Result of @ref computeOrdering(); must outlive this object.
     * @param domains Result of @ref computeDomains(); must outlive this object.
     * @param semantics Whether matches must be induced.
     * @param handler Polled with isRunning() to stop the search on interruption.
     * @param report Receives every complete mapping.
     */
    RIImpl(const SearchGraph &pattern, const SearchGraph &target,
           const std::vector<index> &patternNodeLabels, const std::vector<index> &targetNodeLabels,
           const Ordering &ordering, const Domains &domains,
           SubgraphIsomorphism::Semantics semantics, Aux::SignalHandler &handler,
           MatchReporter report);

    /**
     * @return the empty mapping, with one entry per position in the order.
     */
    State rootState() const;

    /**
     * Enumerates all matches recursively. Used by @ref RI.
     */
    void run();

    /**
     * Expands @a state by one level. Used by @ref ParallelRI. A state at full depth is reported as
     * a match instead.
     *
     * @param state The state to expand. Its `nextCandidate` is advanced past the tried candidates.
     * @param children Receives the children; appended to, not cleared.
     * @return false if the search must stop because the cap on the number of matches is reached.
     */
    bool expand(State &state, std::vector<State> &children);

private:
    /**
     * Maps the pattern node at `state.depth` to every consistent candidate and recurses.
     *
     * @return false if the search must stop because it was interrupted or the cap on the number
     * of matches is reached.
     */
    bool recurse(State &state);

    /**
     * Appends the candidates for position `state.depth` to @a out.
     *
     * If the position has a parent, the candidates are the out- or in-neighbours of the parent's
     * image, depending on the direction of the pattern edge. They are filtered by edge label and,
     * if the domain earns its keep, intersected with the RI-DS domain. Without a parent, the
     * candidates are the RI-DS domain or all target nodes. Target nodes that are already mapped
     * are not filtered out here; @ref ruleEdgesToPrefix() rejects them.
     */
    void candidatesFor(const State &state, std::vector<node> &out) const;

    /**
     * @return true if mapping the pattern node at `state.depth` to @a tv passes the degree check,
     * @ref ruleLabels(), @ref ruleEdgesToPrefix() and, under INDUCED semantics,
     * @ref ruleNonEdgesToPrefix().
     */
    bool consistent(const State &state, node tv) const;

    /**
     * @return true if @a tv is not mapped yet and every pattern arc between the pattern node at
     * `state.depth` and a mapped pattern node corresponds to a target arc with a compatible edge
     * label.
     */
    bool ruleEdgesToPrefix(const State &state, node tv) const;

    /**
     * @return true if no pattern non-edge between the pattern node at `state.depth` and a mapped
     * pattern node corresponds to a target edge. Only called under INDUCED semantics.
     */
    bool ruleNonEdgesToPrefix(const State &state, node tv) const;

    /**
     * @return true if the node labels of @a pu and @a tv are compatible. @ref none is a wildcard
     * on either side.
     */
    bool ruleLabels(node pu, node tv) const;

    /**
     * Reports the complete mapping of @a state. `state.mapping` is indexed by position in the
     * order, so it is permuted into `matchBuffer`, which is indexed by pattern node.
     *
     * @return false if the search must stop.
     */
    bool reportMapping(const State &state);

    const SearchGraph *patternGraph;
    const SearchGraph *targetGraph;

    const std::vector<index> *patternNodeLabels;
    const std::vector<index> *targetNodeLabels;
    bool nodeLabelled;

    const Ordering *ordering;
    /// Empty under plain RI.
    const Domains *domains;

    SubgraphIsomorphism::Semantics semantics;

    /// Only the non-throwing isRunning() may be called, since ParallelRI runs the search inside an
    /// OpenMP region.
    Aux::SignalHandler *handler;

    MatchReporter report;

    /// Reused buffer for reported matches, indexed by pattern node.
    SubgraphIsomorphism::Match matchBuffer;

    /// Candidate lists of all depths of recurse(), stacked in one buffer.
    mutable std::vector<node> candidateBuffer;
};

/**
 * The preprocessing results an RI search needs, as built by @ref prepareRISearch().
 */
struct RISearchSetup {
    /// Snapshot of the pattern, with the adjacency matrix.
    SearchGraph patternGraph;

    /// Snapshot of the target, without the adjacency matrix.
    SearchGraph targetGraph;

    /// Empty under plain RI.
    RIImpl::Domains domains;

    RIImpl::Ordering ordering;
};

/**
 * Builds the snapshots, the RI-DS domains and the matching order for @ref RI and
 * @ref ParallelRI.
 *
 * @param algorithmName Name of the calling algorithm, used in the error message.
 * @throws std::runtime_error if a graph has parallel edges with different edge labels.
 */
RISearchSetup prepareRISearch(const Graph &pattern, const Graph &target,
                              const std::vector<index> &patternNodeLabels,
                              const std::vector<index> &targetNodeLabels,
                              const std::vector<index> &patternEdgeLabels,
                              const std::vector<index> &targetEdgeLabels, RI::Variant variant,
                              const std::string &algorithmName);

} // namespace IsomorphismDetails
} // namespace NetworKit

#endif // NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
