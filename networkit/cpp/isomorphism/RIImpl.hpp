#ifndef NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_

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
 * The RI search, shared by @ref RI and @ref ParallelRI. @ref run() enumerates all matches
 * recursively for RI. @ref expand() produces the children of one @ref State for the workers of
 * ParallelRI, which can steal unexplored states from each other.
 */
class RIImpl {

public:
    /// `order[i]` is the pattern node mapped at depth i. `parent[i]` is the position of an earlier
    /// neighbour of `order[i]`, or @ref none if position i starts a new connected component.
    struct Ordering {
        std::vector<node> order;
        std::vector<index> parent;
    };

    /// RI-DS candidate domains, indexed by pattern node. Empty under plain RI.
    struct Domains {
        /// Sorted target nodes that each pattern node can be mapped to.
        std::vector<std::vector<node>> ofPatternNode;

        /// Whether intersecting candidates with the domain pays off. See
        /// MinSweepYieldForSliceIntersection in RIImpl.cpp.
        std::vector<bool> earnsItsKeep;

        bool anyEmpty = false;
    };

    /// Computes the RI-DS domains from degrees and labels, refines them with one pass over the
    /// pattern arcs, and applies forward checking. Returns empty domains under plain RI.
    static Domains computeDomains(const SearchGraph &pattern, const SearchGraph &target,
                                  const std::vector<index> &patternNodeLabels,
                                  const std::vector<index> &targetNodeLabels, RI::Variant variant);

    /// A node of the search tree. `mapping[i]` is the image of `order[i]` for positions below
    /// `depth`, and `nextCandidate` is where @ref expand() resumes.
    struct State {
        std::vector<node> mapping;
        count depth = 0;
        index nextCandidate = 0;
    };

    /// Computes the matching order of Bonnici et al. Under RI-DS, single-element domains come
    /// first, and a full tie goes to the smaller domain.
    static Ordering computeOrdering(const SearchGraph &pattern, const Domains &domains);

    static bool patternCannotFit(const SearchGraph &pattern, const SearchGraph &target);

    /// @a ordering and @a domains must outlive this object.
    RIImpl(const SearchGraph &pattern, const SearchGraph &target,
           const std::vector<index> &patternNodeLabels, const std::vector<index> &targetNodeLabels,
           const Ordering &ordering, const Domains &domains,
           SubgraphIsomorphism::Semantics semantics, Aux::SignalHandler &handler,
           MatchReporter report);

    State rootState() const;

    void run();

    /// Appends the children of @a state to @a children, or reports @a state if it is complete.
    /// Returns false once the cap on the number of matches is reached.
    bool expand(State &state, std::vector<State> &children);

private:
    bool recurse(State &state);

    /// Appends the candidates for position `state.depth` to @a out. It does not filter out mapped
    /// target nodes; ruleEdgesToPrefix() rejects them.
    void candidatesFor(const State &state, std::vector<node> &out) const;

    bool consistent(const State &state, node tv) const;

    bool ruleEdgesToPrefix(const State &state, node tv) const;

    bool ruleNonEdgesToPrefix(const State &state, node tv) const;

    bool ruleLabels(node pu, node tv) const;

    bool reportMapping(const State &state);

    const SearchGraph *patternGraph;
    const SearchGraph *targetGraph;

    const std::vector<index> *patternNodeLabels;
    const std::vector<index> *targetNodeLabels;
    bool nodeLabelled;

    const Ordering *ordering;
    const Domains *domains;

    SubgraphIsomorphism::Semantics semantics;

    /// Only isRunning() may be called, since ParallelRI runs the search inside an OpenMP region.
    Aux::SignalHandler *handler;

    MatchReporter report;

    SubgraphIsomorphism::Match matchBuffer;

    /// Candidate lists of all depths of recurse(), stacked in one buffer.
    mutable std::vector<node> candidateBuffer;
};

struct RISearchSetup {
    SearchGraph patternGraph;
    SearchGraph targetGraph;
    RIImpl::Domains domains;
    RIImpl::Ordering ordering;
};

/// Builds the snapshots, the RI-DS domains and the matching order. Throws `std::runtime_error` if a
/// graph has parallel edges with different edge labels.
RISearchSetup prepareRISearch(const Graph &pattern, const Graph &target,
                              const std::vector<index> &patternNodeLabels,
                              const std::vector<index> &targetNodeLabels,
                              const std::vector<index> &patternEdgeLabels,
                              const std::vector<index> &targetEdgeLabels, RI::Variant variant,
                              const std::string &algorithmName);

} // namespace IsomorphismDetails
} // namespace NetworKit

#endif // NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
