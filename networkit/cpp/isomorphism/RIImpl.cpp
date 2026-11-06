#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SparseVector.hpp>

#include "RIImpl.hpp"

namespace NetworKit::IsomorphismDetails {

namespace {

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
    const count inArcs = pattern.isDirected() ? pattern.inDegree(pu) : 0;
    out.reserve(out.size() + pattern.outDegree(pu) + inArcs);

    const auto collect = [&](bool outgoing) -> void {
        const Slice slice = sliceOf(pattern, pu, outgoing);
        for (index i = 0; i < slice.nodes.size(); ++i)
            out.push_back(
                {slice.nodes[i], slice.labels.empty() ? none : slice.labels[i], outgoing});
    };

    collect(true);
    if (pattern.isDirected())
        collect(false);
}

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

/// Computes the matching order of Bonnici et al. Under RI-DS, single-element domains come first,
/// and a full tie goes to the smaller domain.
RIOrdering computeOrdering(const SearchGraph &pattern, const RIDomains &domains) {
    // The score of an unordered node is the triple (V_vis, V_neig, V_unv) of Bonnici et al.: its
    // arcs into the ordered prefix, the prefix nodes two hops away through an unordered node, and
    // its unordered neighbors adjacent to no prefix node. Where Figure 2 of the paper disagrees
    // with its text, this follows the text. The best score is reset in every iteration, and ties
    // go to the smallest node id.

    const count z = pattern.upperNodeIdBound();
    const count total = pattern.numberOfNodes();

    std::vector<node> order;
    order.reserve(total);

    std::vector<bool> inOrder(z, false);

    std::vector<count> visCount(z, 0);

    SparseVector<bool> reached(z, false);

    const auto tripleFor = [&](node u) -> std::array<count, 3> {
        std::array<count, 3> score{visCount[u], 0, 0};

        pattern.forEachDistinctNeighbor(u, [&](node v) -> void {
            if (inOrder[v])
                return;

            if (visCount[v] == 0)
                ++score[2];

            pattern.forEachDistinctNeighbor(v, [&](node w) -> void {
                if (!inOrder[w] || reached.indexIsUsed(w))
                    return;
                reached.insert(w, true);
                ++score[1];
            });
        });
        reached.reset();

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

    const auto eligible = [&](node u) -> bool {
        return pattern.hasNode(u) && !inOrder[u]
               && (remainingSingletons == 0 || domainSize(u) == 1);
    };

    while (order.size() < total) {
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

        inOrder[best] = true;
        order.push_back(best);
        pattern.forEachIncidentArc(best, [&](node v) -> void { ++visCount[v]; });
    }

    return orderingFrom(pattern, std::move(order));
}

/// Computes the RI-DS domains from degrees and labels, refines them with one pass over the pattern
/// arcs, and applies forward checking. Returns no domains if one of them becomes empty, in which
/// case there are no matches.
std::optional<RIDomains> computeDomains(const SearchGraph &pattern, const SearchGraph &target,
                                        std::span<const index> patternNodeLabels,
                                        std::span<const index> targetNodeLabels,
                                        double domainFilterThreshold) {
    RIDomains result;
    const count z = pattern.upperNodeIdBound();
    result.ofPatternNode.assign(z, {});
    result.filtersCandidates.assign(z, false);

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

        std::erase_if(domain, [&](node tv) -> bool {
            return std::ranges::any_of(constraints, [&](const ArcConstraint &arc) -> bool {
                return !intersectsDomain(sliceOf(target, tv, arc.outgoing), arc.label,
                                         result.ofPatternNode[arc.pj]);
            });
        });
    }

    if (!forwardCheckSingletons(pattern, result.ofPatternNode))
        return std::nullopt;

    // What the build removed does not count, since consistent() rejects those candidates anyway. A
    // single-node domain always pays off, since the walk over the slice then stops at that node.
    for (node pu = 0; pu < z; ++pu) {
        if (!pattern.hasNode(pu))
            continue;

        const count built = builtSize[pu];
        const count kept = result.ofPatternNode[pu].size();
        const count pruned = built - kept;
        const bool prunedEnough =
            built != 0
            && static_cast<double>(pruned) >= domainFilterThreshold * static_cast<double>(built);
        result.filtersCandidates[pu] = kept == 1 || prunedEnough;
    }

    return result;
}

/// Whether the target is large enough in nodes and degrees to host the pattern. True for an empty
/// pattern.
bool sizesAndDegreesFit(const SearchGraph &pattern, const SearchGraph &target) {
    if (pattern.numberOfNodes() > target.numberOfNodes())
        return false;

    if (pattern.maxOutDegree() > target.maxOutDegree())
        return false;

    return !pattern.isDirected() || pattern.maxInDegree() <= target.maxInDegree();
}

} // namespace

RIOrdering orderingFrom(const SearchGraph &pattern, std::vector<node> order) {
    std::vector<index> positionOf(pattern.upperNodeIdBound(), none);
    for (index i = 0; i < order.size(); ++i)
        positionOf[order[i]] = i;

    RIOrdering result{.order = std::move(order),
                      .parent = {},
                      .checkedNeighborsFirst = {0},
                      .checkedNeighbors = {}};
    const count total = result.order.size();
    result.parent.reserve(total);
    result.checkedNeighborsFirst.reserve(total + 1);

    std::vector<index> &checked = result.checkedNeighbors;
    for (index i = 0; i < total; ++i) {
        const node pu = result.order[i];

        const index begin = checked.size();
        pattern.forEachDistinctNeighbor(pu, [&](node v) -> void {
            if (positionOf[v] < i)
                checked.push_back(positionOf[v]);
        });
        const auto earlier = checked.begin() + static_cast<std::ptrdiff_t>(begin);
        std::sort(earlier, checked.end());

        // A node without earlier neighbors starts a new component and has no parent.
        index parentPos = none;
        if (earlier != checked.end()) {
            parentPos = *earlier;
            const node pp = result.order[parentPos];
            if (!(pattern.isDirected() && pattern.hasEdge(pp, pu) && pattern.hasEdge(pu, pp)))
                checked.erase(earlier);
        }

        result.parent.push_back(parentPos);
        result.checkedNeighborsFirst.push_back(checked.size());
    }

    return result;
}

RISearchSetup prepareRISearch(const Graph &pattern, const Graph &target,
                              std::span<const index> patternNodeLabels,
                              std::span<const index> targetNodeLabels,
                              std::span<const index> patternEdgeLabels,
                              std::span<const index> targetEdgeLabels, RI::Variant variant,
                              double domainFilterThreshold, std::string_view algorithmName) {
    RISearchSetup setup{
        .patternGraph = SearchGraph(pattern, /* buildMatrix = */ true, patternEdgeLabels),
        .targetGraph = SearchGraph(target, /* buildMatrix = */ false, targetEdgeLabels),
        .domains = {},
        .ordering = {}};

    refuseLostLabels(setup.patternGraph, setup.targetGraph, algorithmName);

    // An impossible instance skips the domains and the order.
    setup.canMatch = sizesAndDegreesFit(setup.patternGraph, setup.targetGraph);
    if (!setup.canMatch)
        return setup;

    if (variant == RI::Variant::RI_DS) {
        std::optional<RIDomains> domains =
            computeDomains(setup.patternGraph, setup.targetGraph, patternNodeLabels,
                           targetNodeLabels, domainFilterThreshold);
        setup.canMatch = domains.has_value();
        if (!setup.canMatch)
            return setup;
        setup.domains = std::move(*domains);
    }

    setup.ordering = computeOrdering(setup.patternGraph, setup.domains);
    return setup;
}

} // namespace NetworKit::IsomorphismDetails
