#include <algorithm>
#include <array>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SparseVector.hpp>

#include "RIImpl.hpp"

namespace NetworKit::IsomorphismDetails {

namespace {

/// Ignores direction and visits a directed mutual pair twice, once per arc.
template <typename Callback>
void forEachIncidentArc(const SearchGraph &g, node u, Callback fn) {
    for (const node v : g.outNeighbors(u))
        fn(v);

    // The in-slices of an undirected snapshot are its out-slices.
    if (!g.isDirected())
        return;

    for (const node v : g.inNeighbors(u))
        fn(v);
}

template <typename Callback>
void forEachDistinctNeighbor(const SearchGraph &g, node u, Callback fn) {
    for (const node v : g.outNeighbors(u))
        fn(v);

    if (!g.isDirected())
        return;

    for (const node v : g.inNeighbors(u))
        if (!g.hasEdge(u, v))
            fn(v);
}

bool intersectsDomain(const Slice &slice, index patternLabel, std::span<const node> domain) {
    auto candidate = domain.begin();
    for (index i = 0; i < slice.nodes.size(); ++i) {
        const node tv = slice.nodes[i];
        candidate = std::lower_bound(candidate, domain.end(), tv);

        if (candidate == domain.end())
            return false;

        if (*candidate != tv)
            continue;

        if (slice.labels.empty() || labelsCompatible(patternLabel, slice.labels[i]))
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
    const auto collect = [&](bool outgoing) {
        const Slice slice = sliceOf(pattern, pu, outgoing);
        for (index i = 0; i < slice.nodes.size(); ++i)
            out.push_back(
                {slice.nodes[i], slice.labels.empty() ? none : slice.labels[i], outgoing});
    };

    collect(true);
    if (pattern.isDirected())
        collect(false);
}

/**
 * Minimum fraction of a domain that the refinement and forward checking must remove before the
 * search intersects target neighborhoods with the domain. Below it, the degree and label checks
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
            const auto found = std::ranges::lower_bound(other, fixed);
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

RIOrdering computeOrdering(const SearchGraph &pattern, const RIDomains &domains) {
    // The score of an unordered node is the triple (V_vis, V_neig, V_unv) of Bonnici et al.: its
    // arcs into the ordered prefix, the prefix nodes two hops away through an unordered node, and
    // its unordered neighbors adjacent to no prefix node. Where Figure 2 of the paper disagrees
    // with its text, this follows the text. The best score is reset in every iteration, and ties
    // go to the smallest node id.

    const count z = pattern.upperNodeIdBound();
    const count total = pattern.numberOfNodes();

    RIOrdering result;
    result.order.reserve(total);
    result.parent.reserve(total);

    std::vector<bool> inOrder(z, false);

    std::vector<count> visCount(z, 0);

    // The position of the first ordered neighbor of each node, which becomes its parent.
    std::vector<index> firstOrderedNeighbor(z, none);

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
        const index position = result.order.size();
        inOrder[best] = true;
        result.order.push_back(best);
        result.parent.push_back(firstOrderedNeighbor[best]);
        forEachIncidentArc(pattern, best, [&](node v) {
            ++visCount[v];
            if (firstOrderedNeighbor[v] == none)
                firstOrderedNeighbor[v] = position;
        });
    }

    return result;
}

bool cannotMatch(const SearchGraph &pattern, const SearchGraph &target, const RIDomains &domains) {
    if (domains.anyEmpty || pattern.numberOfNodes() > target.numberOfNodes())
        return true;

    if (pattern.maxOutDegree() > target.maxOutDegree())
        return true;

    return pattern.isDirected() && pattern.maxInDegree() > target.maxInDegree();
}

RIDomains computeDomains(const SearchGraph &pattern, const SearchGraph &target,
                         std::span<const index> patternNodeLabels,
                         std::span<const index> targetNodeLabels, RI::Variant variant) {
    RIDomains result;
    if (variant != RI::Variant::RI_DS)
        return result;

    const count z = pattern.upperNodeIdBound();
    result.ofPatternNode.assign(z, {});
    result.earnsItsKeep.assign(z, false);

    // Walking the target ids in ascending order keeps every domain sorted.
    for (node pu = 0; pu < z; ++pu) {
        if (!pattern.hasNode(pu))
            continue;

        std::vector<node> &domain = result.ofPatternNode[pu];
        for (node tv = 0; tv < target.upperNodeIdBound(); ++tv)
            if (target.hasNode(tv)
                && degreesAndLabelsFit(pattern, target, patternNodeLabels, targetNodeLabels, pu,
                                       tv))
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
            for (const ArcConstraint &arc : constraints)
                if (!intersectsDomain(sliceOf(target, tv, arc.outgoing), arc.label,
                                      result.ofPatternNode[arc.pj]))
                    return false;

            return true;
        };

        std::erase_if(domain, [&](node tv) { return !keep(tv); });
    }

    result.anyEmpty = !forwardCheckSingletons(pattern, result.ofPatternNode);

    // What the build removed does not count, since consistent() rejects those candidates anyway. A
    // single-node domain always pays off, since the walk over the slice then stops at that node.
    for (node pu = 0; pu < z; ++pu) {
        if (!pattern.hasNode(pu))
            continue;

        const count built = builtSize[pu];
        const count kept = result.ofPatternNode[pu].size();
        const count pruned = built - kept;
        result.earnsItsKeep[pu] =
            kept == 1
            || (built != 0
                && static_cast<double>(pruned)
                       >= MinSweepYieldForSliceIntersection * static_cast<double>(built));
    }

    return result;
}

RISearchSetup prepareRISearch(const Graph &pattern, const Graph &target,
                              std::span<const index> patternNodeLabels,
                              std::span<const index> targetNodeLabels,
                              std::span<const index> patternEdgeLabels,
                              std::span<const index> targetEdgeLabels, RI::Variant variant,
                              std::string_view algorithmName) {
    SearchGraph patternGraph(pattern, /* buildMatrix = */ true, patternEdgeLabels);
    SearchGraph targetGraph(target, /* buildMatrix = */ false, targetEdgeLabels);

    refuseLostLabels(patternGraph, targetGraph, algorithmName);

    RIDomains domains =
        computeDomains(patternGraph, targetGraph, patternNodeLabels, targetNodeLabels, variant);
    RIOrdering ordering = computeOrdering(patternGraph, domains);

    return RISearchSetup{std::move(patternGraph), std::move(targetGraph), std::move(domains),
                         std::move(ordering)};
}

} // namespace NetworKit::IsomorphismDetails
