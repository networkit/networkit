#ifndef NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_

#include <algorithm>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/isomorphism/RI.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

#include "MatchReporter.hpp"
#include "SearchGraph.hpp"

namespace NetworKit::IsomorphismDetails {

/// `order[i]` is the pattern node mapped at depth i. `parent[i]` is the position of an earlier
/// neighbor of `order[i]`, or @ref none if position i starts a new connected component.
struct RIOrdering {
    std::vector<node> order;
    std::vector<index> parent;
};

/// RI-DS candidate domains, indexed by pattern node. Empty under plain RI.
struct RIDomains {
    /// Sorted target nodes that each pattern node can be mapped to.
    std::vector<std::vector<node>> ofPatternNode;

    /// Whether intersecting candidates with the domain pays off. See
    /// MinSweepYieldForSliceIntersection in RIImpl.cpp.
    std::vector<bool> earnsItsKeep;

    bool anyEmpty = false;
};

/// A node of the search tree. `mapping[i]` is the image of `order[i]` for positions below
/// `depth`.
struct RIState {
    std::vector<node> mapping;
    count depth = 0;
};

/// Whether the sizes, the degrees or an empty RI-DS domain rule out every match. False for an
/// empty pattern, whose one match is the empty mapping.
bool cannotMatch(const SearchGraph &pattern, const SearchGraph &target, const RIDomains &domains);

/// The out- or in-neighbors of one node, with one arc label per entry.
struct Slice {
    std::span<const node> nodes;
    /// Empty if the graph has no edge labels.
    std::span<const index> labels;
};

// Without `inline`, GCC 15 keeps the call and returns the slice through memory, which slows down
// the search by a few percent.
inline Slice sliceOf(const SearchGraph &g, node u, bool outgoing) {
    if (outgoing)
        return {g.outNeighbors(u), g.outLabels(u)};
    return {g.inNeighbors(u), g.inLabels(u)};
}

/// The checks of a candidate that do not depend on the partial mapping.
inline bool degreesAndLabelsFit(const SearchGraph &pattern, const SearchGraph &target,
                                std::span<const index> patternNodeLabels,
                                std::span<const index> targetNodeLabels, node pu, node tv) {
    if (target.outDegree(tv) < pattern.outDegree(pu))
        return false;

    if (pattern.isDirected() && target.inDegree(tv) < pattern.inDegree(pu))
        return false;

    return nodeLabelsCompatible(patternNodeLabels, targetNodeLabels, pu, tv);
}

/**
 * The RI search behind @ref RI. @ref run() enumerates all matches depth first and passes each one
 * to the reporter, see @ref MatchReporter.
 */
template <MatchReporter Reporter>
class RIImpl {

public:
    /// @a ordering and @a domains must outlive this object.
    RIImpl(const SearchGraph &pattern, const SearchGraph &target,
           std::span<const index> patternNodeLabels, std::span<const index> targetNodeLabels,
           const RIOrdering &ordering, const RIDomains &domains,
           SubgraphIsomorphism::Semantics semantics, Aux::SignalHandler &handler, Reporter report)
        : patternGraph(&pattern), targetGraph(&target), patternNodeLabels(patternNodeLabels),
          targetNodeLabels(targetNodeLabels), ordering(&ordering), domains(&domains),
          semantics(semantics), handler(&handler), report(std::move(report)),
          matchBuffer(pattern.upperNodeIdBound(), none) {}

    RIState rootState() const { return RIState{std::vector<node>(ordering->order.size(), none)}; }

    void run() {
        if (cannotMatch(*patternGraph, *targetGraph, *domains))
            return;

        search();
    }

private:
    /// The depth-first search. It keeps one level of candidates per open depth on the heap, so a
    /// large pattern cannot overflow the call stack. It returns early on an interrupt or once the
    /// cap on the number of matches is reached.
    void search() {
        if (!handler->isRunning())
            return;

        const count matchDepth = ordering->order.size();
        RIState state = rootState();

        // A pattern without nodes has one match, the empty mapping.
        if (matchDepth == 0) {
            reportMapping(state);
            return;
        }

        // The candidates of one open depth are `candidates[begin, end)`. Those from `next` on are
        // still untried.
        struct Level {
            index begin;
            index next;
            index end;
        };

        std::vector<node> candidates;
        std::vector<Level> levels;
        levels.reserve(matchDepth);

        // levels.back() holds the candidates for position state.depth.
        const auto openLevel = [&] {
            const index begin = candidates.size();
            candidatesFor(state, candidates);
            levels.push_back({begin, begin, candidates.size()});
        };

        openLevel();
        while (!levels.empty()) {
            Level &level = levels.back();
            if (level.next == level.end) {
                // We have tried every candidate at this depth, so we go back to the depth above.
                candidates.resize(level.begin);
                levels.pop_back();
                if (!levels.empty())
                    --state.depth;
                continue;
            }

            // Access by index, since deeper levels append to the buffer and may reallocate it.
            const node tv = candidates[level.next++];
            if (!consistent(state, tv))
                continue;

            // Polled only after a candidate passes, since isRunning() is an out-of-line call.
            if (!handler->isRunning())
                return;

            // Nothing reads a mapping entry at or above the depth, so backtracking needs no reset.
            state.mapping[state.depth] = tv;
            if (state.depth + 1 == matchDepth) {
                if (!reportMapping(state))
                    return;
            } else {
                ++state.depth;
                openLevel();
            }
        }
    }

    /// Appends the candidates for position `state.depth` to @a out. It does not filter out mapped
    /// target nodes; ruleEdgesToPrefix() rejects them.
    void candidatesFor(const RIState &state, std::vector<node> &out) const {
        const index pos = state.depth;
        const node pu = ordering->order[pos];
        const index parentPos = ordering->parent[pos];

        const std::vector<node> *builtDomain =
            domains->ofPatternNode.empty() ? nullptr : &domains->ofPatternNode[pu];

        if (parentPos == none) {
            // Without a parent, the domain replaces the scan over all target nodes, whatever its
            // yield.
            if (builtDomain != nullptr) {
                out.insert(out.end(), builtDomain->begin(), builtDomain->end());
                return;
            }

            for (node tv = 0; tv < targetGraph->upperNodeIdBound(); ++tv)
                if (targetGraph->hasNode(tv))
                    out.push_back(tv);
            return;
        }

        const std::vector<node> *domain =
            builtDomain != nullptr && domains->earnsItsKeep[pu] ? builtDomain : nullptr;

        const node pp = ordering->order[parentPos];
        const node parentImage = state.mapping[parentPos];

        // For a mutual pair, walk the shorter slice, since ruleEdgesToPrefix() checks both
        // directions.
        const bool forward = patternGraph->hasEdge(pp, pu);
        const bool backward = patternGraph->isDirected() && patternGraph->hasEdge(pu, pp);
        const bool useOut =
            forward
            && (!backward
                || targetGraph->outDegree(parentImage) <= targetGraph->inDegree(parentImage));
        const Slice slice = sliceOf(*targetGraph, parentImage, useOut);
        const index patternLabel =
            useOut ? patternGraph->edgeLabel(pp, pu) : patternGraph->edgeLabel(pu, pp);

        const bool filterLabels = !slice.labels.empty() && patternLabel != none;

        auto inDomain = domain == nullptr ? std::vector<node>::const_iterator{} : domain->begin();
        for (index i = 0; i < slice.nodes.size(); ++i) {
            const node tv = slice.nodes[i];

            if (filterLabels && !labelsCompatible(patternLabel, slice.labels[i]))
                continue;

            if (domain != nullptr) {
                inDomain = std::lower_bound(inDomain, domain->end(), tv);
                if (inDomain == domain->end())
                    return;
                if (*inDomain != tv)
                    continue;
            }

            out.push_back(tv);
        }
    }

    bool consistent(const RIState &state, node tv) const {
        const node pu = ordering->order[state.depth];

        if (!degreesAndLabelsFit(*patternGraph, *targetGraph, patternNodeLabels, targetNodeLabels,
                                 pu, tv))
            return false;

        if (!ruleEdgesToPrefix(state, tv))
            return false;

        return semantics != SubgraphIsomorphism::Semantics::INDUCED
               || ruleNonEdgesToPrefix(state, tv);
    }

    bool ruleEdgesToPrefix(const RIState &state, node tv) const {
        const node pu = ordering->order[state.depth];
        const bool edgeLabelled = patternGraph->hasEdgeLabels();

        for (index i = 0; i < state.depth; ++i) {
            const node ti = state.mapping[i];

            // Injectivity.
            if (ti == tv)
                return false;

            const node pi = ordering->order[i];

            if (patternGraph->hasEdge(pi, pu)) {
                if (!targetGraph->hasEdge(ti, tv))
                    return false;
                if (edgeLabelled
                    && !labelsCompatible(patternGraph->edgeLabel(pi, pu),
                                         targetGraph->edgeLabel(ti, tv)))
                    return false;
            }

            // The reverse arc of a directed pattern has its own label.
            if (patternGraph->isDirected() && patternGraph->hasEdge(pu, pi)) {
                if (!targetGraph->hasEdge(tv, ti))
                    return false;
                if (edgeLabelled
                    && !labelsCompatible(patternGraph->edgeLabel(pu, pi),
                                         targetGraph->edgeLabel(tv, ti)))
                    return false;
            }
        }

        return true;
    }

    bool ruleNonEdgesToPrefix(const RIState &state, node tv) const {
        const node pu = ordering->order[state.depth];

        for (index i = 0; i < state.depth; ++i) {
            const node pi = ordering->order[i];
            const node ti = state.mapping[i];

            if (!patternGraph->hasEdge(pi, pu) && targetGraph->hasEdge(ti, tv))
                return false;

            if (patternGraph->isDirected() && !patternGraph->hasEdge(pu, pi)
                && targetGraph->hasEdge(tv, ti))
                return false;
        }

        return true;
    }

    bool reportMapping(const RIState &state) {
        for (index i = 0; i < ordering->order.size(); ++i)
            matchBuffer[ordering->order[i]] = state.mapping[i];

        // Every complete mapping writes the same entries, so matchBuffer needs no reset.
        return report(matchBuffer);
    }

    const SearchGraph *patternGraph;
    const SearchGraph *targetGraph;

    std::span<const index> patternNodeLabels;
    std::span<const index> targetNodeLabels;

    const RIOrdering *ordering;
    const RIDomains *domains;

    SubgraphIsomorphism::Semantics semantics;

    /// The search only polls isRunning() and returns. The caller throws for the interrupt.
    Aux::SignalHandler *handler;

    Reporter report;

    std::vector<node> matchBuffer;
};

struct RISearchSetup {
    SearchGraph patternGraph;
    SearchGraph targetGraph;
    RIDomains domains;
    RIOrdering ordering;
};

/// Builds the snapshots, the RI-DS domains and the matching order. Throws `std::runtime_error` if a
/// graph has parallel edges with different edge labels.
RISearchSetup prepareRISearch(const Graph &pattern, const Graph &target,
                              std::span<const index> patternNodeLabels,
                              std::span<const index> targetNodeLabels,
                              std::span<const index> patternEdgeLabels,
                              std::span<const index> targetEdgeLabels, RI::Variant variant,
                              std::string_view algorithmName);

} // namespace NetworKit::IsomorphismDetails

#endif // NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
