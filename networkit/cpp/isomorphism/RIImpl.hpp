#ifndef NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_

#include <algorithm>
#include <optional>
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

/// `order[i]` is the pattern node mapped at depth i. `parent[i]` is the position of the first
/// earlier neighbor of `order[i]`, or @ref none if position i starts a new connected component.
/// @ref orderingFrom() fills the other fields from these two.
struct RIOrdering {
    std::vector<node> order;
    std::vector<index> parent;

    /// See @ref checkedNeighborsOf().
    std::vector<index> checkedNeighborsFirst;
    std::vector<index> checkedNeighbors;

    /// The positions before @a i whose pattern nodes share an arc with `order[i]`, in ascending
    /// order. The parent is left out, since the candidates for position i come from the target
    /// slice along its arc, unless it forms a directed mutual pair with `order[i]`. The candidates
    /// then come from only one of the two arcs.
    std::span<const index> checkedNeighborsOf(index i) const {
        return std::span<const index>(checkedNeighbors)
            .subspan(checkedNeighborsFirst[i],
                     checkedNeighborsFirst[i + 1] - checkedNeighborsFirst[i]);
    }
};

/// Completes @a order to an RIOrdering. Every node of @a pattern must occur in @a order once.
RIOrdering orderingFrom(const SearchGraph &pattern, std::vector<node> order);

/// RI-DS candidate domains, indexed by pattern node. Empty under plain RI.
struct RIDomains {
    /// Sorted target nodes that each pattern node can be mapped to.
    std::vector<std::vector<node>> ofPatternNode;

    /// Whether candidatesFor() filters the candidates of a pattern node through its domain. See
    /// the parameter `domainFilterThreshold` of @ref RI.
    std::vector<bool> filtersCandidates;
};

/// A node of the search tree. `mapping[i]` is the image of `order[i]` for positions below
/// `depth`.
struct RIState {
    std::vector<node> mapping;
    count depth = 0;
};

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

/// What @ref prepareRISearch() computes before the search.
struct RISearchSetup {
    SearchGraph patternGraph;
    SearchGraph targetGraph;
    RIDomains domains;
    RIOrdering ordering;

    /// False if the sizes, the degrees or an empty RI-DS domain rule out every match. The domains
    /// and the order are then not computed. True for an empty pattern, whose one match is the
    /// empty mapping.
    bool canMatch = true;
};

/**
 * The RI search, shared by @ref RI and @ref ParallelRI. @ref run() enumerates all matches depth
 * first for RI. @ref expand() produces the children of one @ref RIState for the workers of
 * ParallelRI, which can steal unexplored states from each other. Both pass each match to the
 * reporter, see @ref MatchReporter.
 */
template <MatchReporter Reporter>
class RIImpl {

public:
    /// @a setup must outlive this object. The buffers are sized here rather than in run(), because
    /// the workers of ParallelRI only call expand().
    RIImpl(const RISearchSetup &setup, std::span<const index> patternNodeLabels,
           std::span<const index> targetNodeLabels, SubgraphIsomorphism::Semantics semantics,
           Aux::SignalHandler &handler, Reporter report)
        : patternGraph(&setup.patternGraph), targetGraph(&setup.targetGraph),
          patternNodeLabels(patternNodeLabels), targetNodeLabels(targetNodeLabels),
          ordering(&setup.ordering), domains(&setup.domains), canMatch(setup.canMatch),
          semantics(semantics), handler(&handler), report(std::move(report)),
          matchBuffer(setup.patternGraph.upperNodeIdBound(), none),
          mapped(setup.targetGraph.upperNodeIdBound(), false) {}

    RIState rootState() const { return RIState{std::vector<node>(ordering->order.size(), none)}; }

    void run() {
        if (canMatch)
            search();
    }

    /// Appends the children of @a state to @a children, or reports @a state if it is complete.
    /// Children that complete the mapping are reported right away instead, except for those of the
    /// root. Like search(), it polls for an interrupt first and before every match. @a state must
    /// descend from rootState(). Returns false on an interrupt or once the cap on the number of
    /// matches is reached. Unlike run(), it does not check RISearchSetup::canMatch.
    bool expand(const RIState &state, std::vector<RIState> &children) {
        if (!handler->isRunning())
            return false;

        const count matchDepth = ordering->order.size();
        if (state.depth == matchDepth)
            return reportMapping(state);

        candidateBuffer.clear();
        candidatesFor(state, candidateBuffer);

        // consistent() reads the marks of the mapped target nodes, which a state does not carry.
        // The guard clears them on every way out, also if the reporter throws.
        markPrefix(state, true);
        const struct ClearMarks {
            RIImpl &impl;
            const RIState &state;
            ~ClearMarks() { impl.markPrefix(state, false); }
        } clearMarks{*this, state};

        // A child that completes the mapping is reported here, which saves a state per match. The
        // children of the root stay states, so that ParallelRI deals them out even for a one-node
        // pattern.
        const bool reportChildren = state.depth != 0 && state.depth + 1 == matchDepth;
        if (reportChildren)
            for (index i = 0; i < state.depth; ++i)
                matchBuffer[ordering->order[i]] = state.mapping[i];

        for (const node tv : candidateBuffer) {
            if (!consistent(state, tv))
                continue;

            if (reportChildren) {
                // A slow callback must not hold up an interrupt for a whole level of candidates.
                if (!handler->isRunning())
                    return false;

                matchBuffer[ordering->order[state.depth]] = tv;
                if (!report(matchBuffer))
                    return false;
                continue;
            }

            RIState child = state;
            child.mapping[state.depth] = tv;
            child.depth = state.depth + 1;
            children.push_back(std::move(child));
        }

        return true;
    }

private:
    /// Sets the marks of the target nodes that @a state maps to @a mark.
    void markPrefix(const RIState &state, bool mark) {
        for (index i = 0; i < state.depth; ++i)
            mapped[state.mapping[i]] = mark;
    }

    /// The depth-first search. It keeps one level of candidates per open depth on the heap, so a
    /// large pattern cannot overflow the call stack. It returns early on an interrupt or once the
    /// cap on the number of matches is reached.
    void search() {
        if (!handler->isRunning())
            return;

        // An earlier run may have stopped with target nodes still marked.
        mapped.assign(mapped.size(), false);

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
        const auto openLevel = [&]() -> void {
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
                if (!levels.empty()) {
                    --state.depth;
                    mapped[state.mapping[state.depth]] = false;
                }
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
                mapped[tv] = true;
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

        if (parentPos == none) {
            // Without a parent, the domain replaces the scan over all target nodes, whatever its
            // yield.
            if (!domains->ofPatternNode.empty()) {
                const std::vector<node> &domain = domains->ofPatternNode[pu];
                out.insert(out.end(), domain.begin(), domain.end());
                return;
            }

            for (node tv = 0; tv < targetGraph->upperNodeIdBound(); ++tv)
                if (targetGraph->hasNode(tv))
                    out.push_back(tv);
            return;
        }

        const bool useDomain = !domains->ofPatternNode.empty() && domains->filtersCandidates[pu];
        std::span<const node> domain;
        if (useDomain)
            domain = domains->ofPatternNode[pu];

        const node pp = ordering->order[parentPos];
        const node parentImage = state.mapping[parentPos];

        // Every candidate has the arc of the walked slice, so ruleEdgesToPrefix() does not check
        // it again. For a mutual pair, walk the shorter slice, since ruleEdgesToPrefix() checks
        // both directions.
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

        auto inDomain = domain.begin();
        for (index i = 0; i < slice.nodes.size(); ++i) {
            const node tv = slice.nodes[i];

            if (filterLabels && !labelsCompatible(patternLabel, slice.labels[i]))
                continue;

            if (useDomain) {
                inDomain = std::lower_bound(inDomain, domain.end(), tv);
                if (inDomain == domain.end())
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

    /// Checks injectivity and the arcs to the positions of RIOrdering::checkedNeighborsOf().
    bool ruleEdgesToPrefix(const RIState &state, node tv) const {
        if (mapped[tv])
            return false;

        const node pu = ordering->order[state.depth];

        for (const index i : ordering->checkedNeighborsOf(state.depth)) {
            const node pi = ordering->order[i];
            const node ti = state.mapping[i];

            if (patternGraph->hasEdge(pi, pu) && !arcFits(pi, pu, ti, tv))
                return false;

            // The reverse arc of a directed pattern has its own label.
            if (patternGraph->isDirected() && patternGraph->hasEdge(pu, pi)
                && !arcFits(pu, pi, tv, ti))
                return false;
        }

        return true;
    }

    /// Whether the target has the arc @a tFrom -> @a tTo with a label that fits the pattern arc
    /// @a pFrom -> @a pTo.
    bool arcFits(node pFrom, node pTo, node tFrom, node tTo) const {
        if (!patternGraph->hasEdgeLabels())
            return targetGraph->hasEdge(tFrom, tTo);

        // One search finds both the arc and its label.
        const std::optional<index> targetLabel = targetGraph->findArc(tFrom, tTo);
        return targetLabel.has_value()
               && labelsCompatible(patternGraph->edgeLabel(pFrom, pTo), *targetLabel);
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

    /// See @ref RISearchSetup::canMatch.
    bool canMatch;

    SubgraphIsomorphism::Semantics semantics;

    /// Only isRunning() may be called, since ParallelRI calls expand() inside an OpenMP region. The
    /// caller throws for the interrupt.
    Aux::SignalHandler *handler;

    Reporter report;

    std::vector<node> matchBuffer;

    /// Marks the target nodes of the mapped prefix, for the injectivity check.
    std::vector<bool> mapped;

    /// The candidates of the state that expand() works on. A member, so that the allocation
    /// carries over from one call to the next.
    std::vector<node> candidateBuffer;
};

/// Builds the snapshots, decides whether a match is possible, and if so, computes the RI-DS
/// domains and the matching order. See @ref RI for @a domainFilterThreshold. Throws
/// `std::runtime_error` if a graph has parallel edges with different edge labels.
RISearchSetup prepareRISearch(const Graph &pattern, const Graph &target,
                              std::span<const index> patternNodeLabels,
                              std::span<const index> targetNodeLabels,
                              std::span<const index> patternEdgeLabels,
                              std::span<const index> targetEdgeLabels, RI::Variant variant,
                              double domainFilterThreshold, std::string_view algorithmName);

} // namespace NetworKit::IsomorphismDetails

#endif // NETWORKIT_CPP_ISOMORPHISM_RI_IMPL_HPP_
