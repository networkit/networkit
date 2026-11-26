#ifndef NETWORKIT_CPP_ISOMORPHISM_TEST_SUBGRAPH_ISOMORPHISM_TEST_UTILS_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_TEST_SUBGRAPH_ISOMORPHISM_TEST_UTILS_HPP_

// Shared test helpers for the isomorphism module: a brute-force reference matcher, a corpus of
// test cases, and assertions that compare an algorithm against the reference.

#include <algorithm>
#include <initializer_list>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit::IsomorphismTest {

using Match = SubgraphIsomorphism::Match;
using Semantics = SubgraphIsomorphism::Semantics;

// The reference matcher works on Graph rather than on SearchGraph, so that it cannot share a bug
// with the code under test.

/// The label @ref none matches any label.
inline bool labelsCompatible(index patternLabel, index targetLabel) {
    return patternLabel == none || targetLabel == none || patternLabel == targetLabel;
}

inline bool nodeLabelsCompatible(std::span<const index> patternNodeLabels,
                                 std::span<const index> targetNodeLabels, node pu, node tv) {
    return patternNodeLabels.empty()
           || labelsCompatible(patternNodeLabels[pu], targetNodeLabels[tv]);
}

/// The labels of all edges, keyed by ordered node pair. Parallel edges give a pair several labels.
using LabelsByPair = std::map<std::pair<node, node>, std::vector<index>>;

inline LabelsByPair edgeLabelsByPair(const Graph &G, std::span<const index> edgeLabels) {
    LabelsByPair byPair;
    G.forEdges([&](node u, node v, edgeweight, edgeid eid) -> void {
        const index label = eid < edgeLabels.size() ? edgeLabels[eid] : none;
        byPair[{u, v}].push_back(label);
        if (!G.isDirected() && u != v)
            byPair[{v, u}].push_back(label);
    });
    return byPair;
}

inline const std::vector<index> &labelsOfPair(const LabelsByPair &byPair, node u, node v) {
    static const std::vector<index> nothing;
    const auto found = byPair.find({u, v});
    return found == byPair.end() ? nothing : found->second;
}

/// With parallel edges, some pattern label must be compatible with some target label.
inline bool edgeLabelsCompatible(std::span<const index> patternEdgeLabels,
                                 std::span<const index> targetEdgeLabels) {
    return std::ranges::any_of(patternEdgeLabels, [&](index patternLabel) -> bool {
        return std::ranges::any_of(targetEdgeLabels, [&](index targetLabel) -> bool {
            return labelsCompatible(patternLabel, targetLabel);
        });
    });
}

/// Whether a complete mapping is a well-formed match.
inline bool isValidMatch(const Graph &pattern, const Graph &target, Semantics semantics,
                         std::span<const index> patternNodeLabels,
                         std::span<const index> targetNodeLabels, std::span<const node> match,
                         std::span<const index> patternEdgeLabels = {},
                         std::span<const index> targetEdgeLabels = {}) {
    const count pz = pattern.upperNodeIdBound();

    if (match.size() != pz)
        return false;

    const bool edgeLabelled = !patternEdgeLabels.empty();
    LabelsByPair patternEdgeLabelsByPair;
    LabelsByPair targetEdgeLabelsByPair;
    if (edgeLabelled) {
        patternEdgeLabelsByPair = edgeLabelsByPair(pattern, patternEdgeLabels);
        targetEdgeLabelsByPair = edgeLabelsByPair(target, targetEdgeLabels);
    }

    for (node pu = 0; pu < pz; ++pu) {
        if (!pattern.hasNode(pu)) {
            if (match[pu] != none)
                return false;
            continue;
        }
        if (match[pu] == none || match[pu] >= target.upperNodeIdBound()
            || !target.hasNode(match[pu]))
            return false;
        if (!nodeLabelsCompatible(patternNodeLabels, targetNodeLabels, pu, match[pu]))
            return false;
    }

    // Injective: no two pattern nodes share an image.
    std::vector<node> images;
    for (const node u : pattern.nodeRange())
        images.push_back(match[u]);
    std::ranges::sort(images);
    if (std::ranges::adjacent_find(images) != images.end())
        return false;

    // Ordered pairs check both directions of a directed graph.
    for (const node a : pattern.nodeRange()) {
        for (const node b : pattern.nodeRange()) {
            if (a == b)
                continue;

            const bool patternEdge = pattern.hasEdge(a, b);
            const bool targetEdge = target.hasEdge(match[a], match[b]);

            if (patternEdge && !targetEdge)
                return false;
            if (semantics == Semantics::INDUCED && !patternEdge && targetEdge)
                return false;

            // Edge labels only constrain pattern edges.
            if (edgeLabelled && patternEdge
                && !edgeLabelsCompatible(labelsOfPair(patternEdgeLabelsByPair, a, b),
                                         labelsOfPair(targetEdgeLabelsByPair, match[a], match[b])))
                return false;
        }
    }

    return true;
}

/// Checks every injective mapping, which takes O(targetNodes ^ patternNodes) time.
inline std::vector<Match> referenceMatches(const Graph &pattern, const Graph &target,
                                           Semantics semantics,
                                           std::span<const index> patternNodeLabels = {},
                                           std::span<const index> targetNodeLabels = {},
                                           std::span<const index> patternEdgeLabels = {},
                                           std::span<const index> targetEdgeLabels = {}) {
    const std::vector<node> patternNodes(pattern.nodeRange().begin(), pattern.nodeRange().end());

    std::vector<Match> matches;
    Match current(pattern.upperNodeIdBound(), none);
    std::vector<bool> used(target.upperNodeIdBound(), false);

    const auto assign = [&](const auto &self, index next) -> void {
        if (next == patternNodes.size()) {
            if (isValidMatch(pattern, target, semantics, patternNodeLabels, targetNodeLabels,
                             current, patternEdgeLabels, targetEdgeLabels))
                matches.push_back(current);
            return;
        }

        const node pu = patternNodes[next];
        for (const node tv : target.nodeRange()) {
            if (used[tv])
                continue;

            used[tv] = true;
            current[pu] = tv;
            self(self, next + 1);
            current[pu] = none;
            used[tv] = false;
        }
    };
    assign(assign, 0);

    return matches;
}

/// Edge label vectors are indexed by edge id, so their graphs come from labelledGraphOf(). The
/// labels default to empty, so that a case built with designated initializers can leave them out.
struct Case {
    std::string name;
    Graph pattern;
    Graph target;
    Semantics semantics;
    std::vector<index> patternNodeLabels{};
    std::vector<index> targetNodeLabels{};
    std::vector<index> patternEdgeLabels{};
    std::vector<index> targetEdgeLabels{};
};

inline Case unlabelledCase(std::string name, Graph pattern, Graph target, Semantics semantics) {
    return {.name = std::move(name),
            .pattern = std::move(pattern),
            .target = std::move(target),
            .semantics = semantics};
}

inline Graph graphOf(count n, std::initializer_list<std::pair<node, node>> edges,
                     bool directed = false) {
    Graph G(n, false, directed);
    for (const auto &[u, v] : edges)
        G.addEdge(u, v);
    return G;
}

/// The complete graph on @a n nodes. The directed one has both arcs between every pair.
inline Graph completeGraph(count n, bool directed = false) {
    Graph G(n, false, directed);
    for (node u = 0; u < n; ++u)
        for (node v = directed ? 0 : u + 1; v < n; ++v)
            if (u != v)
                G.addEdge(u, v);
    return G;
}

struct LabelledGraph {
    Graph G;
    std::vector<index> edgeLabels;
};

/// Like graphOf(), with each edge written as `{u, v, label}`. Parallel edges may swap labels.
inline LabelledGraph labelledGraphOf(count n,
                                     std::initializer_list<std::tuple<node, node, index>> edges,
                                     bool directed = false) {
    const auto key = [directed](node u, node v) -> std::pair<node, node> {
        return directed || u <= v ? std::pair<node, node>{u, v} : std::pair<node, node>{v, u};
    };

    Graph G(n, false, directed);
    LabelsByPair queued;
    for (const auto &[u, v, label] : edges) {
        G.addEdge(u, v);
        queued[key(u, v)].push_back(label);
    }
    G.indexEdges();

    std::vector<index> edgeLabels(G.upperEdgeIdBound(), none);
    std::map<std::pair<node, node>, index> handedOut;
    G.forEdges([&](node u, node v, edgeweight, edgeid eid) -> void {
        const std::pair<node, node> pair = key(u, v);
        edgeLabels[eid] = queued[pair][handedOut[pair]++];
    });

    return {std::move(G), std::move(edgeLabels)};
}

inline std::vector<Case> standardCases() {
    std::vector<Case> cases;

    const Graph triangle = graphOf(3, {{0, 1}, {1, 2}, {2, 0}});
    const Graph path3 = graphOf(3, {{0, 1}, {1, 2}});
    const Graph k4 = completeGraph(4);
    const Graph k5 = completeGraph(5);

    // A 3-path occurs in a triangle, but not as an induced occurrence.
    cases.push_back(
        unlabelledCase("path3-in-triangle-mono", path3, triangle, Semantics::MONOMORPHISM));
    cases.push_back(
        unlabelledCase("path3-in-triangle-induced", path3, triangle, Semantics::INDUCED));

    cases.push_back(unlabelledCase("triangle-in-k4-mono", triangle, k4, Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("triangle-in-k4-induced", triangle, k4, Semantics::INDUCED));
    cases.push_back(unlabelledCase("k4-in-k5-induced", k4, k5, Semantics::INDUCED));

    cases.push_back(unlabelledCase("pattern-larger-than-target", triangle, graphOf(2, {{0, 1}}),
                                   Semantics::MONOMORPHISM));

    // An isolated pattern node takes its candidates from all target nodes.
    cases.push_back(unlabelledCase("pattern-with-isolated-node", graphOf(3, {{0, 1}}),
                                   graphOf(4, {{0, 1}, {1, 2}}), Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("disconnected-pattern", graphOf(4, {{0, 1}, {2, 3}}),
                                   graphOf(5, {{0, 1}, {1, 2}, {3, 4}}), Semantics::MONOMORPHISM));

    // A removed id has an empty slice, like an isolated node.
    Graph gappedTarget = graphOf(7, {{0, 1}, {1, 3}, {3, 4}, {4, 6}});
    gappedTarget.removeNode(2);
    gappedTarget.removeNode(5);
    cases.push_back(
        unlabelledCase("target-with-removed-ids", path3, gappedTarget, Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("isolated-pattern-node-vs-removed-ids", graphOf(3, {{0, 1}}),
                                   gappedTarget, Semantics::MONOMORPHISM));

    Graph gappedPattern = graphOf(4, {{0, 1}, {1, 3}});
    gappedPattern.removeNode(2);
    cases.push_back(
        unlabelledCase("pattern-with-removed-ids", gappedPattern, k4, Semantics::MONOMORPHISM));

    // Both degenerate patterns have one match, the empty mapping.
    cases.push_back(unlabelledCase("empty-pattern", Graph(0), k4, Semantics::INDUCED));
    Graph allRemoved(3);
    for (node u = 0; u < 3; ++u)
        allRemoved.removeNode(u);
    cases.push_back(
        unlabelledCase("pattern-with-all-nodes-removed", allRemoved, k4, Semantics::INDUCED));
    cases.push_back(
        unlabelledCase("single-node-pattern", Graph(1), gappedTarget, Semantics::MONOMORPHISM));

    Graph messyTarget(4, false, false);
    messyTarget.addEdge(0, 1);
    messyTarget.addEdge(0, 1); // parallel
    messyTarget.addEdge(1, 2);
    messyTarget.addEdge(2, 0);
    messyTarget.addEdge(3, 3); // self-loop
    messyTarget.addEdge(2, 3);
    cases.push_back(unlabelledCase("target-with-multiedges-and-loop", path3, messyTarget,
                                   Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("target-with-multiedges-and-loop-induced", path3, messyTarget,
                                   Semantics::INDUCED));

    // The 2-cycle catches an implementation that checks only one direction.
    const Graph arc = graphOf(2, {{0, 1}}, true);
    const Graph twoCycle = graphOf(2, {{0, 1}, {1, 0}}, true);
    const Graph directedTriangle = graphOf(3, {{0, 1}, {1, 2}, {2, 0}}, true);
    const Graph directedTarget = graphOf(4, {{0, 1}, {1, 0}, {1, 2}, {2, 3}, {3, 1}}, true);

    cases.push_back(
        unlabelledCase("arc-in-directed", arc, directedTarget, Semantics::MONOMORPHISM));
    cases.push_back(
        unlabelledCase("two-cycle-in-directed", twoCycle, directedTarget, Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("directed-triangle", directedTriangle, directedTarget,
                                   Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("directed-triangle-induced", directedTriangle, directedTarget,
                                   Semantics::INDUCED));

    // Node 0 has only in-arcs, so VF2 continues from the in-terminal sets.
    const Graph inArcsFirst = graphOf(3, {{1, 0}, {2, 0}, {2, 1}}, true);
    const Graph completeDigraph = completeGraph(4, true);

    cases.push_back(
        unlabelledCase("in-arcs-first-induced", inArcsFirst, directedTarget, Semantics::INDUCED));
    cases.push_back(unlabelledCase("in-arcs-first-in-complete-digraph", inArcsFirst,
                                   completeDigraph, Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("directed-empty-pattern", Graph(0, false, true), directedTarget,
                                   Semantics::INDUCED));

    // Each star exceeds the target's maximum degree in one direction, but not its size.
    const Graph oneEdge = graphOf(5, {{1, 2}});
    const Graph oneArc = graphOf(5, {{1, 2}}, true);
    cases.push_back(unlabelledCase("star-exceeds-max-degree", graphOf(4, {{0, 1}, {0, 2}, {0, 3}}),
                                   oneEdge, Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("in-star-exceeds-max-in-degree",
                                   graphOf(4, {{1, 0}, {2, 0}, {3, 0}}, true), oneArc,
                                   Semantics::MONOMORPHISM));
    cases.push_back(unlabelledCase("out-star-exceeds-max-out-degree",
                                   graphOf(4, {{0, 1}, {0, 2}, {0, 3}}, true), oneArc,
                                   Semantics::MONOMORPHISM));

    cases.push_back({.name = "labelled-path3-in-k4",
                     .pattern = path3,
                     .target = k4,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternNodeLabels = {7, 8, 7},
                     .targetNodeLabels = {7, 8, 7, 9}});
    cases.push_back({.name = "labelled-wildcard-on-pattern",
                     .pattern = path3,
                     .target = k4,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternNodeLabels = {none, 8, none},
                     .targetNodeLabels = {7, 8, 7, 9}});
    cases.push_back({.name = "labelled-wildcard-on-target",
                     .pattern = path3,
                     .target = k4,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternNodeLabels = {7, 8, 7},
                     .targetNodeLabels = {none, 8, 7, none}});
    cases.push_back({.name = "labelled-no-match",
                     .pattern = triangle,
                     .target = k4,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternNodeLabels = {1, 1, 1},
                     .targetNodeLabels = {1, 1, 2, 2}});

    // Most edge-label cases use MONOMORPHISM, since the induced rule constrains non-edges, which
    // carry no labels.
    const LabelledGraph edgeLabelledPattern = labelledGraphOf(3, {{0, 1, 1}, {1, 2, 2}});
    const LabelledGraph edgeLabelledTarget =
        labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}, {0, 3, 3}});

    cases.push_back({.name = "edge-labelled-match",
                     .pattern = edgeLabelledPattern.G,
                     .target = edgeLabelledTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = edgeLabelledPattern.edgeLabels,
                     .targetEdgeLabels = edgeLabelledTarget.edgeLabels});

    // An implementation that ignores edge labels reports matches here.
    const LabelledGraph nearMissTarget =
        labelledGraphOf(4, {{0, 1, 1}, {1, 2, 5}, {2, 3, 1}, {0, 3, 3}});
    cases.push_back({.name = "edge-labelled-near-miss",
                     .pattern = edgeLabelledPattern.G,
                     .target = nearMissTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = edgeLabelledPattern.edgeLabels,
                     .targetEdgeLabels = nearMissTarget.edgeLabels});

    const LabelledGraph wildcardPattern = labelledGraphOf(3, {{0, 1, none}, {1, 2, 2}});
    cases.push_back({.name = "edge-labelled-wildcard-on-pattern",
                     .pattern = wildcardPattern.G,
                     .target = edgeLabelledTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = wildcardPattern.edgeLabels,
                     .targetEdgeLabels = edgeLabelledTarget.edgeLabels});

    const LabelledGraph wildcardTarget =
        labelledGraphOf(4, {{0, 1, none}, {1, 2, 2}, {2, 3, 1}, {0, 3, none}});
    cases.push_back({.name = "edge-labelled-wildcard-on-target",
                     .pattern = edgeLabelledPattern.G,
                     .target = wildcardTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = edgeLabelledPattern.edgeLabels,
                     .targetEdgeLabels = wildcardTarget.edgeLabels});

    // The two arcs of a directed mutual pair carry independent labels.
    const LabelledGraph mutualPattern = labelledGraphOf(2, {{0, 1, 7}}, true);
    const LabelledGraph mutualTarget = labelledGraphOf(3, {{0, 1, 7}, {1, 0, 8}, {1, 2, 8}}, true);
    cases.push_back({.name = "directed-mutual-pair-different-labels",
                     .pattern = mutualPattern.G,
                     .target = mutualTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = mutualPattern.edgeLabels,
                     .targetEdgeLabels = mutualTarget.edgeLabels});

    // A search may walk one arc of a mutual pair and only look up the other. From target nodes 0
    // and 2 such a search walks different arcs, and each time the other arc has the wrong label.
    const LabelledGraph twoLabelPair = labelledGraphOf(2, {{0, 1, 1}, {1, 0, 2}}, true);
    const LabelledGraph twoLabelTarget = labelledGraphOf(
        6, {{0, 1, 1}, {1, 0, 5}, {3, 0, 2}, {2, 3, 9}, {3, 2, 2}, {2, 1, 4}, {4, 5, 1}, {5, 4, 2}},
        true);
    cases.push_back({.name = "directed-mutual-pair-labels-checked-both-ways",
                     .pattern = twoLabelPair.G,
                     .target = twoLabelTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = twoLabelPair.edgeLabels,
                     .targetEdgeLabels = twoLabelTarget.edgeLabels});

    // Parallel edges with equal labels collapse without loss.
    const LabelledGraph parallelAgreeing =
        labelledGraphOf(4, {{0, 1, 1}, {0, 1, 1}, {1, 2, 2}, {2, 3, 1}});
    cases.push_back({.name = "edge-labelled-parallel-edges-same-label",
                     .pattern = edgeLabelledPattern.G,
                     .target = parallelAgreeing.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = edgeLabelledPattern.edgeLabels,
                     .targetEdgeLabels = parallelAgreeing.edgeLabels});

    const LabelledGraph parallelAgreeingPattern =
        labelledGraphOf(3, {{0, 1, 1}, {0, 1, 1}, {1, 2, 2}});
    cases.push_back({.name = "edge-labelled-parallel-pattern-edges-same-label",
                     .pattern = parallelAgreeingPattern.G,
                     .target = edgeLabelledTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = parallelAgreeingPattern.edgeLabels,
                     .targetEdgeLabels = edgeLabelledTarget.edgeLabels});

    // The snapshot drops target self-loops, so their labels cannot disagree.
    const LabelledGraph labelledLoops =
        labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}, {2, 2, 7}, {2, 2, 8}});
    cases.push_back({.name = "edge-labelled-target-self-loops",
                     .pattern = edgeLabelledPattern.G,
                     .target = labelledLoops.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternEdgeLabels = edgeLabelledPattern.edgeLabels,
                     .targetEdgeLabels = labelledLoops.edgeLabels});

    // The node labels rule out one of the two edge-labelled matches.
    cases.push_back({.name = "node-and-edge-labelled",
                     .pattern = edgeLabelledPattern.G,
                     .target = edgeLabelledTarget.G,
                     .semantics = Semantics::MONOMORPHISM,
                     .patternNodeLabels = {7, none, 8},
                     .targetNodeLabels = {7, 9, 8, 7},
                     .patternEdgeLabels = edgeLabelledPattern.edgeLabels,
                     .targetEdgeLabels = edgeLabelledTarget.edgeLabels});

    // The labelled square occurs in the target, but the diagonal makes every occurrence
    // non-induced.
    const LabelledGraph labelledSquare =
        labelledGraphOf(4, {{0, 1, 0}, {1, 2, 42}, {2, 3, 0}, {3, 0, 42}});
    const LabelledGraph squareWithDiagonal =
        labelledGraphOf(4, {{0, 1, 42}, {1, 2, 0}, {2, 3, 42}, {3, 0, 0}, {0, 2, 42}});
    for (const Semantics semantics : {Semantics::MONOMORPHISM, Semantics::INDUCED})
        cases.push_back({.name = semantics == Semantics::INDUCED ? "edge-labelled-square-induced"
                                                                 : "edge-labelled-square-mono",
                         .pattern = labelledSquare.G,
                         .target = squareWithDiagonal.G,
                         .semantics = semantics,
                         .patternEdgeLabels = labelledSquare.edgeLabels,
                         .targetEdgeLabels = squareWithDiagonal.edgeLabels});

    // Seeded random graphs reach deeper backtracking than the small cases above. The directed ones
    // catch a search that maps two pattern nodes to one target node through the in-terminal sets,
    // or that skips the induced check of the in-arcs.
    Aux::Random::setSeed(42, false);
    const auto randomGraph = [](count n, bool directed) -> Graph {
        Graph G(n, false, directed);
        for (node u = 0; u < n; ++u)
            for (node v = directed ? 0 : u + 1; v < n; ++v)
                if (u != v && Aux::Random::probability() < 0.5)
                    G.addEdge(u, v);
        return G;
    };
    for (index i = 0; i < 40; ++i) {
        const bool directed = i % 2 == 1;
        cases.push_back(unlabelledCase("random-" + std::to_string(i), randomGraph(4, directed),
                                       randomGraph(8, directed),
                                       i % 4 < 2 ? Semantics::INDUCED : Semantics::MONOMORPHISM));
    }

    return cases;
}

/// The reference matcher behind the interface, following the sequential run() protocol.
class ReferenceSubgraphIsomorphism final : public SubgraphIsomorphism {

public:
    ReferenceSubgraphIsomorphism(const Graph &pattern, const Graph &target,
                                 Semantics semantics = Semantics::INDUCED, count maxMatches = 0)
        : SubgraphIsomorphism(pattern, target, semantics, maxMatches) {}

    void run() override {
        Aux::SignalHandler handler;

        prepareRun();

        for (const Match &match :
             referenceMatches(*pattern, *target, semantics, patternNodeLabels, targetNodeLabels,
                              patternEdgeLabels, targetEdgeLabels)) {
            handler.assureRunning();
            if (!reportMatch(match))
                break;
        }

        finishRun();
    }
};

// The assertions take a factory rather than a type, so that an algorithm with more constructor
// parameters can use them.

inline void applyLabels(SubgraphIsomorphism &algo, const Case &testCase) {
    if (!testCase.patternNodeLabels.empty())
        algo.setNodeLabels(testCase.patternNodeLabels, testCase.targetNodeLabels);
    if (!testCase.patternEdgeLabels.empty())
        algo.setEdgeLabels(testCase.patternEdgeLabels, testCase.targetEdgeLabels);
}

/// Runs every standard case with stored matches, with a callback and with counting only, and
/// compares each run with the reference.
template <typename Construct>
void expectMatchesReference(Construct construct) {
    for (const Case &testCase : standardCases()) {
        SCOPED_TRACE("case: " + testCase.name);

        std::vector<Match> expected = referenceMatches(
            testCase.pattern, testCase.target, testCase.semantics, testCase.patternNodeLabels,
            testCase.targetNodeLabels, testCase.patternEdgeLabels, testCase.targetEdgeLabels);
        // UnorderedElementsAreArray would compare every pair of matches, which is slow for the
        // random cases with over a thousand matches. The checks below compare sorted lists instead,
        // and WhenSorted() sorts only the matches of the algorithm.
        std::ranges::sort(expected);

        const auto build = [&]() -> std::unique_ptr<SubgraphIsomorphism> {
            std::unique_ptr<SubgraphIsomorphism> algo =
                construct(testCase.pattern, testCase.target, testCase.semantics, count{0});
            applyLabels(*algo, testCase);
            return algo;
        };

        {
            std::unique_ptr<SubgraphIsomorphism> algo = build();
            algo->run();

            EXPECT_THAT(algo->getMatches(),
                        testing::WhenSorted(testing::ElementsAreArray(expected)));
            EXPECT_EQ(algo->numberOfMatches(), expected.size());
            EXPECT_EQ(algo->hasMatch(), !expected.empty());
        }

        // The serial callback needs no lock, since it is never called concurrently.
        {
            std::vector<Match> collected;
            std::unique_ptr<SubgraphIsomorphism> algo = build();
            algo->setCallback([&](std::span<const node> match) -> void {
                collected.emplace_back(match.begin(), match.end());
            });
            algo->run();

            EXPECT_THAT(collected, testing::WhenSorted(testing::ElementsAreArray(expected)))
                << "serial callback";
            EXPECT_EQ(algo->numberOfMatches(), expected.size()) << "serial callback";
        }

        {
            std::unique_ptr<SubgraphIsomorphism> algo = build();
            algo->setStoreMatches(false);
            algo->run();

            EXPECT_EQ(algo->numberOfMatches(), expected.size()) << "counting only";
            EXPECT_THROW(algo->getMatches(), std::runtime_error) << "counting only";
        }
    }
}

template <typename Construct>
void expectRespectsMatchCap(Construct construct) {
    // 24 matches: every ordering of 3 of K4's 4 nodes.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}, {2, 0}});
    const Graph target = completeGraph(4);
    const count total = 24;

    for (count cap : {count{1}, count{5}, count{23}, total, total + 10}) {
        std::unique_ptr<SubgraphIsomorphism> algo =
            construct(pattern, target, Semantics::MONOMORPHISM, cap);
        algo->run();

        const count expected = std::min(cap, total);
        EXPECT_EQ(algo->numberOfMatches(), expected) << "cap: " << cap;
        EXPECT_THAT(algo->getMatches(), testing::SizeIs(expected)) << "cap: " << cap;
        EXPECT_TRUE(algo->hasMatch()) << "cap: " << cap;

        for (const Match &match : algo->getMatches()) {
            EXPECT_TRUE(isValidMatch(pattern, target, Semantics::MONOMORPHISM, {}, {}, match))
                << "cap: " << cap << " produced a malformed match";
        }
    }
}

/// A callback receives exactly as many matches as the cap allows, and an exception from it stops
/// run() and reaches the caller. A parallel search may pass a callback a few matches beyond the
/// cap, so this suits only the sequential algorithms.
template <typename Construct>
void expectCallbackStopsTheSearch(Construct construct) {
    // 24 matches: every ordering of 3 of K4's 4 nodes.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}, {2, 0}});
    const Graph target = completeGraph(4);
    const count total = 24;

    for (count cap : {count{1}, count{5}, total, total + 10}) {
        std::vector<Match> collected;
        std::unique_ptr<SubgraphIsomorphism> algo =
            construct(pattern, target, Semantics::MONOMORPHISM, cap);
        algo->setCallback([&](std::span<const node> match) -> void {
            collected.emplace_back(match.begin(), match.end());
        });
        algo->run();

        EXPECT_THAT(collected, testing::SizeIs(std::min(cap, total))) << "cap: " << cap;
        EXPECT_EQ(algo->numberOfMatches(), collected.size()) << "cap: " << cap;
        for (const Match &match : collected) {
            EXPECT_TRUE(isValidMatch(pattern, target, Semantics::MONOMORPHISM, {}, {}, match))
                << "cap: " << cap << " produced a malformed match";
        }
    }

    std::unique_ptr<SubgraphIsomorphism> algo =
        construct(pattern, target, Semantics::MONOMORPHISM, count{0});
    count delivered = 0;
    algo->setCallback([&](std::span<const node>) -> void {
        if (++delivered == 3)
            throw std::runtime_error("the callback gave up");
    });

    EXPECT_THAT([&]() -> void { algo->run(); },
                testing::ThrowsMessage<std::runtime_error>(testing::HasSubstr("gave up")));
    EXPECT_EQ(delivered, 3u) << "the search went on after the callback threw";
    EXPECT_FALSE(algo->hasFinished());
    EXPECT_THROW(algo->numberOfMatches(), std::runtime_error);
}

} // namespace NetworKit::IsomorphismTest

#endif // NETWORKIT_CPP_ISOMORPHISM_TEST_SUBGRAPH_ISOMORPHISM_TEST_UTILS_HPP_
