#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SparseVector.hpp>

#include "RIImpl.hpp"

namespace NetworKit {
namespace IsomorphismDetails {

namespace {

/// @ref none is a wildcard on either side.
bool nodeLabelsCompatible(const std::vector<index> &patternNodeLabels,
                          const std::vector<index> &targetNodeLabels, node pu, node tv) {
    if (patternNodeLabels.empty())
        return true;

    const index patternLabel = patternNodeLabels[pu];
    const index targetLabel = targetNodeLabels[tv];
    return patternLabel == none || targetLabel == none || patternLabel == targetLabel;
}

bool edgeLabelsCompatible(index patternLabel, index targetLabel) noexcept {
    return patternLabel == none || targetLabel == none || patternLabel == targetLabel;
}

/// Ignores direction and visits a directed mutual pair twice, once per arc.
template <typename Callback>
void forEachIncidentArc(const SearchGraph &g, node u, Callback fn) {
    for (const node *it = g.outBegin(u); it != g.outEnd(u); ++it)
        fn(*it);

    // The in-slices of an undirected snapshot are its out-slices.
    if (!g.isDirected())
        return;

    for (const node *it = g.inBegin(u); it != g.inEnd(u); ++it)
        fn(*it);
}

template <typename Callback>
void forEachDistinctNeighbor(const SearchGraph &g, node u, Callback fn) {
    for (const node *it = g.outBegin(u); it != g.outEnd(u); ++it)
        fn(*it);

    if (!g.isDirected())
        return;

    for (const node *it = g.inBegin(u); it != g.inEnd(u); ++it)
        if (!g.hasEdge(u, *it))
            fn(*it);
}

bool couldMap(const SearchGraph &pattern, const SearchGraph &target,
              const std::vector<index> &patternNodeLabels,
              const std::vector<index> &targetNodeLabels, node pu, node tv) {
    if (!target.hasNode(tv))
        return false;

    if (target.outDegree(tv) < pattern.outDegree(pu))
        return false;

    if (pattern.isDirected() && target.inDegree(tv) < pattern.inDegree(pu))
        return false;

    return nodeLabelsCompatible(patternNodeLabels, targetNodeLabels, pu, tv);
}

/// @a labels holds one label per slice entry, or nullptr to ignore labels.
bool intersectsDomain(const node *begin, const node *end, const index *labels, index patternLabel,
                      const std::vector<node> &domain) {
    auto candidate = domain.begin();
    for (const node *it = begin; it != end; ++it) {
        candidate = std::lower_bound(candidate, domain.end(), *it);

        if (candidate == domain.end())
            return false;

        if (*candidate != *it)
            continue;

        if (labels == nullptr || edgeLabelsCompatible(patternLabel, labels[it - begin]))
            return true;
    }

    return false;
}

/// One pattern arc at the pattern node whose domain the refinement pass filters.
struct ArcConstraint {
    node pj;
    index label;
    bool outgoing;
};

/// A directed mutual pair yields two constraints, since its two arcs may carry different labels.
void collectArcConstraints(const SearchGraph &pattern, node pu, std::vector<ArcConstraint> &out) {
    const node *outBegin = pattern.outBegin(pu);
    const index *outLabels = pattern.outLabelBegin(pu);
    for (const node *it = outBegin; it != pattern.outEnd(pu); ++it)
        out.push_back({*it, outLabels == nullptr ? none : outLabels[it - outBegin], true});

    if (!pattern.isDirected())
        return;

    const node *inBegin = pattern.inBegin(pu);
    const index *inLabels = pattern.inLabelBegin(pu);
    for (const node *it = inBegin; it != pattern.inEnd(pu); ++it)
        out.push_back({*it, inLabels == nullptr ? none : inLabels[it - inBegin], false});
}

/**
 * Minimum fraction of a domain that the refinement and forward checking must remove before the
 * search intersects target neighbourhoods with the domain. Below it, the degree and label checks
 * in RIImpl::consistent() reject the same candidates more cheaply. On caidaRouterLevel with a
 * labelled triangle pattern, intersecting costs 21% at a yield of 0.56 and breaks even at about
 * 0.8.
 */
constexpr double MinSweepYieldForSliceIntersection = 0.8;

/// Forward checking, Section 4.2.2 of Kimmig, Meyerhenke and Strash. Returns false if a domain
/// becomes empty, in which case there are no matches.
bool forwardCheckSingletons(const SearchGraph &pattern, std::vector<std::vector<node>> &domains) {
    const count z = pattern.upperNodeIdBound();

    std::vector<bool> queued(z, false);
    std::vector<node> pending;

    for (node pu = 0; pu < z; ++pu) {
        if (!pattern.hasNode(pu))
            continue;
        if (domains[pu].empty())
            return false;
        if (domains[pu].size() == 1) {
            queued[pu] = true;
            pending.push_back(pu);
        }
    }

    while (!pending.empty()) {
        const node pu = pending.back();
        pending.pop_back();

        const node fixed = domains[pu].front();

        for (node pw = 0; pw < z; ++pw) {
            if (pw == pu || !pattern.hasNode(pw))
                continue;

            std::vector<node> &other = domains[pw];
            const auto found = std::lower_bound(other.begin(), other.end(), fixed);
            if (found == other.end() || *found != fixed)
                continue;

            other.erase(found);

            // Two pattern nodes need the same target node.
            if (other.empty())
                return false;

            if (other.size() == 1 && !queued[pw]) {
                queued[pw] = true;
                pending.push_back(pw);
            }
        }
    }

    return true;
}

} // namespace

RIImpl::Ordering RIImpl::computeOrdering(const SearchGraph &pattern, const Domains &domains) {
    // The score of an unordered node is the triple (V_vis, V_neig, V_unv) of Bonnici et al.: its
    // arcs into the ordered prefix, the prefix nodes two hops away through an unordered node, and
    // its unordered neighbours adjacent to no prefix node. Where Figure 2 of the paper disagrees
    // with its text, this follows the text. The best score is reset in every iteration, and ties
    // go to the smallest node id.

    const count z = pattern.upperNodeIdBound();
    const count total = pattern.numberOfNodes();

    Ordering result;
    result.order.reserve(total);
    result.parent.reserve(total);

    std::vector<bool> inOrder(z, false);

    std::vector<count> visCount(z, 0);

    SparseVector<bool> reached(z, false);

    const auto tripleFor = [&](node u) {
        std::array<count, 3> score{visCount[u], 0, 0};

        forEachDistinctNeighbor(pattern, u, [&](node v) {
            if (inOrder[v])
                return;

            forEachDistinctNeighbor(pattern, v, [&](node w) {
                if (!inOrder[w] || reached.indexIsUsed(w))
                    return;
                reached.insert(w, true);
                ++score[1];
            });
        });
        reached.reset();

        forEachDistinctNeighbor(pattern, u, [&](node v) {
            if (!inOrder[v] && visCount[v] == 0)
                ++score[2];
        });

        return score;
    };

    // 0 under plain RI, which disables both RI-DS rules below.
    const auto domainSize = [&](node u) -> count {
        return domains.ofPatternNode.empty() ? 0 : domains.ofPatternNode[u].size();
    };

    // Under RI-DS, single-element domains come first (Section 4.1 of Kimmig, Meyerhenke and
    // Strash).
    count remainingSingletons = 0;
    for (node u = 0; u < z; ++u)
        if (pattern.hasNode(u) && domainSize(u) == 1)
            ++remainingSingletons;

    const auto eligible = [&](node u) {
        return pattern.hasNode(u) && !inOrder[u]
               && (remainingSingletons == 0 || domainSize(u) == 1);
    };

    while (result.order.size() < total) {
        // Compute the expensive triple only for the nodes that are best on the first term.
        count bestVis = 0;
        for (node u = 0; u < z; ++u)
            if (eligible(u))
                bestVis = std::max(bestVis, visCount[u]);

        node best = none;
        std::array<count, 3> bestScore{};
        for (node u = 0; u < z; ++u) {
            if (!eligible(u) || visCount[u] != bestVis)
                continue;

            const std::array<count, 3> score = tripleFor(u);

            // Under RI-DS, a tie on the triple goes to the smaller domain first (Section 4.2.1).
            if (best == none || score > bestScore
                || (score == bestScore && domainSize(u) < domainSize(best))) {
                best = u;
                bestScore = score;
            }
        }

        if (remainingSingletons != 0)
            --remainingSingletons;

        // A node without arcs into the order starts a new component and has no parent.
        index parentPos = none;
        if (visCount[best] != 0) {
            for (index j = 0; j < result.order.size(); ++j) {
                const node earlier = result.order[j];
                if (pattern.hasEdge(earlier, best) || pattern.hasEdge(best, earlier)) {
                    parentPos = j;
                    break;
                }
            }
        }

        inOrder[best] = true;
        result.order.push_back(best);
        result.parent.push_back(parentPos);
        forEachIncidentArc(pattern, best, [&](node v) { ++visCount[v]; });
    }

    return result;
}

bool RIImpl::patternCannotFit(const SearchGraph &pattern, const SearchGraph &target) {
    if (pattern.numberOfNodes() > target.numberOfNodes())
        return true;

    if (pattern.maxOutDegree() > target.maxOutDegree())
        return true;

    return pattern.isDirected() && pattern.maxInDegree() > target.maxInDegree();
}

RIImpl::RIImpl(const SearchGraph &pattern, const SearchGraph &target,
               const std::vector<index> &patternNodeLabels,
               const std::vector<index> &targetNodeLabels, const Ordering &ordering,
               const Domains &domains, SubgraphIsomorphism::Semantics semantics,
               Aux::SignalHandler &handler, MatchReporter report)
    : patternGraph(&pattern), targetGraph(&target), patternNodeLabels(&patternNodeLabels),
      targetNodeLabels(&targetNodeLabels), nodeLabelled(!patternNodeLabels.empty()),
      ordering(&ordering), domains(&domains), semantics(semantics), handler(&handler),
      report(std::move(report)) {

    // Sized here rather than in run(), because the workers of ParallelRI only call expand().
    matchBuffer.assign(pattern.upperNodeIdBound(), none);
}

RIImpl::State RIImpl::rootState() const {
    return State{std::vector<node>(ordering->order.size(), none)};
}

void RIImpl::run() {
    // An empty pattern has exactly one match, the empty mapping, so it must not return early.
    if (!ordering->order.empty()) {
        if (patternCannotFit(*patternGraph, *targetGraph))
            return;

        if (domains->anyEmpty)
            return;
    }

    State state = rootState();
    recurse(state);
}

bool RIImpl::recurse(State &state) {
    if (!handler->isRunning())
        return false;

    if (state.depth == ordering->order.size())
        return reportMapping(state);

    const index base = candidateBuffer.size();
    candidatesFor(state, candidateBuffer);
    const index end = candidateBuffer.size();

    bool keepGoing = true;
    for (index k = base; k < end && keepGoing; ++k) {
        // Access by index, since deeper levels append to the buffer and may reallocate it.
        const node tv = candidateBuffer[k];
        if (!consistent(state, tv))
            continue;

        state.mapping[state.depth] = tv;
        ++state.depth;
        keepGoing = recurse(state);
        --state.depth;
        state.mapping[state.depth] = none;
    }

    candidateBuffer.resize(base);
    return keepGoing;
}

bool RIImpl::expand(State &state, std::vector<State> &children) {
    const count full = ordering->order.size();

    if (state.depth == full)
        return reportMapping(state);

    candidateBuffer.clear();
    candidatesFor(state, candidateBuffer);

    for (index k = state.nextCandidate; k < candidateBuffer.size(); ++k) {
        // Advance first, so that a resumed expansion continues after this candidate.
        state.nextCandidate = k + 1;

        const node tv = candidateBuffer[k];
        if (!consistent(state, tv))
            continue;

        State child = state;
        child.mapping[state.depth] = tv;
        child.depth = state.depth + 1;
        child.nextCandidate = 0;
        children.push_back(std::move(child));
    }

    state.nextCandidate = candidateBuffer.size();

    return true;
}

void RIImpl::candidatesFor(const State &state, std::vector<node> &out) const {
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

    // For a mutual pair, walk the shorter slice, since ruleEdgesToPrefix() checks both directions.
    const bool forward = patternGraph->hasEdge(pp, pu);
    const bool backward = patternGraph->isDirected() && patternGraph->hasEdge(pu, pp);
    bool useOut = forward;
    if (forward && backward) {
        useOut = targetGraph->outDegree(parentImage) <= targetGraph->inDegree(parentImage);
    }

    const node *begin =
        useOut ? targetGraph->outBegin(parentImage) : targetGraph->inBegin(parentImage);
    const node *end = useOut ? targetGraph->outEnd(parentImage) : targetGraph->inEnd(parentImage);
    const index *sliceLabels =
        useOut ? targetGraph->outLabelBegin(parentImage) : targetGraph->inLabelBegin(parentImage);
    const index patternLabel =
        useOut ? patternGraph->edgeLabel(pp, pu) : patternGraph->edgeLabel(pu, pp);

    const bool filterLabels = sliceLabels != nullptr && patternLabel != none;

    auto inDomain = domain == nullptr ? std::vector<node>::const_iterator{} : domain->begin();
    for (const node *it = begin; it != end; ++it) {
        const node tv = *it;

        if (filterLabels && !edgeLabelsCompatible(patternLabel, sliceLabels[it - begin]))
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

bool RIImpl::consistent(const State &state, node tv) const {
    const node pu = ordering->order[state.depth];

    if (targetGraph->outDegree(tv) < patternGraph->outDegree(pu))
        return false;

    if (patternGraph->isDirected() && targetGraph->inDegree(tv) < patternGraph->inDegree(pu))
        return false;

    if (!ruleLabels(pu, tv))
        return false;

    if (!ruleEdgesToPrefix(state, tv))
        return false;

    return semantics != SubgraphIsomorphism::Semantics::INDUCED || ruleNonEdgesToPrefix(state, tv);
}

bool RIImpl::ruleEdgesToPrefix(const State &state, node tv) const {
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
                && !edgeLabelsCompatible(patternGraph->edgeLabel(pi, pu),
                                         targetGraph->edgeLabel(ti, tv)))
                return false;
        }

        // The reverse arc of a directed pattern has its own label.
        if (patternGraph->isDirected() && patternGraph->hasEdge(pu, pi)) {
            if (!targetGraph->hasEdge(tv, ti))
                return false;
            if (edgeLabelled
                && !edgeLabelsCompatible(patternGraph->edgeLabel(pu, pi),
                                         targetGraph->edgeLabel(tv, ti)))
                return false;
        }
    }

    return true;
}

bool RIImpl::ruleNonEdgesToPrefix(const State &state, node tv) const {
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

bool RIImpl::ruleLabels(node pu, node tv) const {
    return nodeLabelsCompatible(*patternNodeLabels, *targetNodeLabels, pu, tv);
}

RIImpl::Domains RIImpl::computeDomains(const SearchGraph &pattern, const SearchGraph &target,
                                       const std::vector<index> &patternNodeLabels,
                                       const std::vector<index> &targetNodeLabels,
                                       RI::Variant variant) {
    Domains result;
    if (variant != RI::Variant::RI_DS)
        return result;

    const count z = pattern.upperNodeIdBound();
    if (z == 0)
        return result;

    result.ofPatternNode.assign(z, {});
    result.earnsItsKeep.assign(z, false);

    // Walking the target ids in ascending order keeps every domain sorted.
    for (node pu = 0; pu < z; ++pu) {
        if (!pattern.hasNode(pu))
            continue;

        std::vector<node> &domain = result.ofPatternNode[pu];
        for (node tv = 0; tv < target.upperNodeIdBound(); ++tv)
            if (couldMap(pattern, target, patternNodeLabels, targetNodeLabels, pu, tv))
                domain.push_back(tv);
    }

    // A single refinement pass, since the pass is costly on a large unlabelled target, where a
    // domain holds almost every node. Refining in place is sound, because a removed node cannot be
    // the image of its pattern node in any match.
    std::vector<count> builtSize(z, 0);
    std::vector<ArcConstraint> constraints;

    for (node pu = 0; pu < z; ++pu) {
        if (!pattern.hasNode(pu))
            continue;

        std::vector<node> &domain = result.ofPatternNode[pu];
        builtSize[pu] = domain.size();

        constraints.clear();
        collectArcConstraints(pattern, pu, constraints);

        const auto keep = [&](node tv) {
            for (const ArcConstraint &arc : constraints) {
                const node *begin = arc.outgoing ? target.outBegin(tv) : target.inBegin(tv);
                const node *end = arc.outgoing ? target.outEnd(tv) : target.inEnd(tv);
                const index *labels = arc.label == none ? nullptr
                                                        : (arc.outgoing ? target.outLabelBegin(tv)
                                                                        : target.inLabelBegin(tv));

                if (!intersectsDomain(begin, end, labels, arc.label, result.ofPatternNode[arc.pj]))
                    return false;
            }

            return true;
        };

        domain.erase(
            std::remove_if(domain.begin(), domain.end(), [&](node tv) { return !keep(tv); }),
            domain.end());
    }

    result.anyEmpty = !forwardCheckSingletons(pattern, result.ofPatternNode);

    // What the build removed does not count, since consistent() rejects those candidates anyway.
    for (node pu = 0; pu < z; ++pu) {
        if (!pattern.hasNode(pu))
            continue;

        const count built = builtSize[pu];
        const count pruned = built - result.ofPatternNode[pu].size();
        result.earnsItsKeep[pu] =
            built != 0
            && static_cast<double>(pruned)
                   >= MinSweepYieldForSliceIntersection * static_cast<double>(built);
    }

    return result;
}

bool RIImpl::reportMapping(const State &state) {
    for (index i = 0; i < ordering->order.size(); ++i)
        matchBuffer[ordering->order[i]] = state.mapping[i];

    // Every complete mapping writes the same entries, so matchBuffer needs no reset.
    return report(matchBuffer);
}

RISearchSetup prepareRISearch(const Graph &pattern, const Graph &target,
                              const std::vector<index> &patternNodeLabels,
                              const std::vector<index> &targetNodeLabels,
                              const std::vector<index> &patternEdgeLabels,
                              const std::vector<index> &targetEdgeLabels, RI::Variant variant,
                              const std::string &algorithmName) {
    SearchGraph patternGraph(pattern, /* buildMatrix = */ true, patternEdgeLabels);
    SearchGraph targetGraph(target, /* buildMatrix = */ false, targetEdgeLabels);

    if (patternGraph.collapsedLabelledEdges() || targetGraph.collapsedLabelledEdges())
        throw std::runtime_error(algorithmName
                                 + " does not support parallel edges whose edge labels disagree - "
                                   "see SubgraphIsomorphism::setEdgeLabels()");

    RIImpl::Domains domains = RIImpl::computeDomains(patternGraph, targetGraph, patternNodeLabels,
                                                     targetNodeLabels, variant);
    RIImpl::Ordering ordering = RIImpl::computeOrdering(patternGraph, domains);

    return RISearchSetup{std::move(patternGraph), std::move(targetGraph), std::move(domains),
                         std::move(ordering)};
}

} // namespace IsomorphismDetails
} // namespace NetworKit
