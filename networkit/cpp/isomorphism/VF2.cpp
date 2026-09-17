#include <array>
#include <stdexcept>
#include <vector>

#include <tlx/unused.hpp>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/isomorphism/VF2.hpp>

#include "MatchReporter.hpp"
#include "SearchGraph.hpp"

namespace NetworKit {

namespace {

using IsomorphismDetails::MatchReporter;
using IsomorphismDetails::SearchGraph;

class VF2Impl {

public:
    VF2Impl(const Graph &pattern, const Graph &target, const std::vector<index> &patternNodeLabels,
            const std::vector<index> &targetNodeLabels, const std::vector<index> &patternEdgeLabels,
            const std::vector<index> &targetEdgeLabels, SubgraphIsomorphism::Semantics semantics,
            Aux::SignalHandler &handler, MatchReporter report)
        : patternGraph(pattern, /* buildMatrix = */ true, patternEdgeLabels),
          targetGraph(target, /* buildMatrix = */ false, targetEdgeLabels),
          patternNodeLabels(&patternNodeLabels), targetNodeLabels(&targetNodeLabels),
          nodeLabelled(!patternNodeLabels.empty()), edgeLabelled(!patternEdgeLabels.empty()),
          semantics(semantics), handler(&handler), report(std::move(report)), t1in(0), t1out(0),
          t2in(0), t2out(0) {
        if (patternGraph.collapsedLabelledEdges()) {
            throw std::runtime_error(
                "VF2 does not run if pattern has unequally-labelled collapsed edges.");
        }
        if (targetGraph.collapsedLabelledEdges()) {
            throw std::runtime_error(
                "VF2 does not run if target has unequally-labelled collapsed edges.");
        }
    }

    void run() {

        core1.assign(patternGraph.upperNodeIdBound(), none);
        core2.assign(targetGraph.upperNodeIdBound(), none);
        in1.assign(patternGraph.upperNodeIdBound(), none);
        out1.assign(patternGraph.upperNodeIdBound(), none);
        in2.assign(targetGraph.upperNodeIdBound(), none);
        out2.assign(targetGraph.upperNodeIdBound(), none);
        mapping.resize(patternGraph.upperNodeIdBound(), none);

        if (patternGraph.numberOfNodes() == 0) {
            reportMapping();
            return;
        }

        if (patternGraph.numberOfNodes() > targetGraph.numberOfNodes()
            || patternGraph.maxInDegree() > targetGraph.maxInDegree()
            || patternGraph.maxOutDegree() > targetGraph.maxOutDegree()) {
            return;
        }

        match(0);
    }

private:
    /// One level of the depth-first search. Returns false if the whole search must stop.
    bool match(count depth) {

        if (depth == patternGraph.numberOfNodes()) {
            return reportMapping();
        }

        index cursor = 0;
        node pu = none;
        node tv = none;
        bool continueSearch;

        while (nextCandidatePair(cursor, pu, tv)) {
            handler->assureRunning();
            if (feasible(pu, tv)) {
                const std::array<count, 8> restoreTerminalSets = addPair(pu, tv);
                continueSearch = match(depth + 1);
                removePair(pu, tv, restoreTerminalSets);
                if (!continueSearch) {
                    return false;
                }
            }
        }

        return true;
    }

    /// Produces the next candidate pair at this depth. @a cursor records where the previous call
    /// stopped. Returns false once the candidates are exhausted.
    bool nextCandidatePair(index &cursor, node &pu, node &tv) const {

        if (t1out != 0 && t2out != 0) {

            pu = smallestMember(membersOut1, out1);
            if (pu == none) {
                return false;
            }
            for (index i = cursor; i < membersOut2.size(); ++i) {
                node v = membersOut2[i];
                if (out2[v] != none) {
                    tv = v;
                    cursor = i + 1;
                    return true;
                }
            }

            return false;

        } else if (t1in != 0 && t2in != 0) {

            pu = smallestMember(membersIn1, in1);
            if (pu == none) {
                return false;
            }
            for (index i = cursor; i < membersIn2.size(); ++i) {
                node v = membersIn2[i];
                if (in2[v] != none) {
                    tv = v;
                    cursor = i + 1;
                    return true;
                }
            }

            return false;

        } else {

            for (node u = 0; u < core1.size(); ++u) {
                if (patternGraph.hasNode(u) && core1[u] == none) {
                    pu = u;
                    break;
                }
            }
            for (node v = cursor; v < core2.size(); ++v) {
                if (targetGraph.hasNode(v) && core2[v] == none) {
                    tv = v;
                    cursor = v + 1;
                    return true;
                }
            }
        }

        return false;
    }

    bool feasible(node pu, node tv) const {

        return ruleSuccessors(pu, tv) && rulePredecessors(pu, tv) && ruleTerminalCounts(pu, tv)
               && ruleNewCounts(pu, tv) && ruleLabels(pu, tv);
    }

    /// Every mapped out-neighbour of @a pu needs a corresponding out-edge of @a tv with a
    /// compatible edge label. Under INDUCED semantics, the same holds in reverse.
    bool ruleSuccessors(node pu, node tv) const {

        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it) {
            node u = *it;
            if (core1[u] != none) {
                if (!targetGraph.hasEdge(tv, core1[u])) {
                    return false;
                }
                if (edgeLabelled) {
                    if (!(patternGraph.edgeLabel(pu, u) == targetGraph.edgeLabel(tv, core1[u])
                          || patternGraph.edgeLabel(pu, u) == none
                          || targetGraph.edgeLabel(tv, core1[u]) == none)) {
                        return false;
                    }
                }
            }
        }

        if (semantics == SubgraphIsomorphism::Semantics::INDUCED) {
            for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
                node v = *it;
                if (core2[v] != none) {
                    if (!patternGraph.hasEdge(pu, core2[v])) {
                        return false;
                    }
                    if (edgeLabelled) {
                        if (!(patternGraph.edgeLabel(pu, core2[v]) == targetGraph.edgeLabel(tv, v)
                              || patternGraph.edgeLabel(pu, core2[v]) == none
                              || targetGraph.edgeLabel(tv, v) == none)) {
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    /// The mirror image of ruleSuccessors() for in-edges.
    bool rulePredecessors(node pu, node tv) const {

        if (!patternGraph.isDirected()) {
            return true;
        }

        for (auto it = patternGraph.inBegin(pu); it != patternGraph.inEnd(pu); ++it) {
            node u = *it;
            if (core1[u] != none) {
                if (!targetGraph.hasEdge(core1[u], tv)) {
                    return false;
                }
                if (edgeLabelled) {
                    if (!(patternGraph.edgeLabel(u, pu) == targetGraph.edgeLabel(core1[u], tv)
                          || patternGraph.edgeLabel(u, pu) == none
                          || targetGraph.edgeLabel(core1[u], tv) == none)) {
                        return false;
                    }
                }
            }
        }

        if (semantics == SubgraphIsomorphism::Semantics::INDUCED) {
            for (auto it = targetGraph.inBegin(tv); it != targetGraph.inEnd(tv); ++it) {
                node v = *it;
                if (core2[v] != none) {
                    if (!patternGraph.hasEdge(core2[v], pu)) {
                        return false;
                    }
                    if (edgeLabelled) {
                        if (!(patternGraph.edgeLabel(core2[v], pu) == targetGraph.edgeLabel(v, tv)
                              || patternGraph.edgeLabel(core2[v], pu) == none
                              || targetGraph.edgeLabel(v, tv) == none)) {
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    /// One-step look-ahead: @a pu must not have more unmapped neighbours in the terminal sets than
    /// @a tv has.
    bool ruleTerminalCounts(node pu, node tv) const {

        count in1Neighbors = 0;
        count out1Neighbors = 0;
        count in2Neighbors = 0;
        count out2Neighbors = 0;

        for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
            node v = *it;
            if (core2[v] == none && out2[v] != none) {
                out2Neighbors++;
            }
        }
        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it) {
            node u = *it;
            if (core1[u] == none && out1[u] != none) {
                if (++out1Neighbors > out2Neighbors) {
                    return false;
                }
            }
        }

        if (patternGraph.isDirected()) {
            for (auto it = targetGraph.inBegin(tv); it != targetGraph.inEnd(tv); ++it) {
                node v = *it;
                if (core2[v] == none && in2[v] != none) {
                    in2Neighbors++;
                }
            }
            for (auto it = patternGraph.inBegin(pu); it != patternGraph.inEnd(pu); ++it) {
                node u = *it;
                if (core1[u] == none && in1[u] != none) {
                    if (++in1Neighbors > in2Neighbors) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    /// Two-step look-ahead: the same count for unmapped neighbours outside all terminal sets.
    bool ruleNewCounts(node pu, node tv) const {

        // MONOMORPHISM allows extra target edges, so this rule does not apply.
        if (semantics == SubgraphIsomorphism::Semantics::MONOMORPHISM) {
            return true;
        }

        count in1Neighbors = 0;
        count out1Neighbors = 0;
        count in2Neighbors = 0;
        count out2Neighbors = 0;

        for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
            node v = *it;
            if (core2[v] == none && in2[v] == none && out2[v] == none) {
                out2Neighbors++;
            }
        }
        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it) {
            node u = *it;
            if (core1[u] == none && in1[u] == none && out1[u] == none) {
                if (++out1Neighbors > out2Neighbors) {
                    return false;
                }
            }
        }

        if (patternGraph.isDirected()) {
            for (auto it = targetGraph.inBegin(tv); it != targetGraph.inEnd(tv); ++it) {
                node v = *it;
                if (core2[v] == none && in2[v] == none && out2[v] == none) {
                    in2Neighbors++;
                }
            }
            for (auto it = patternGraph.inBegin(pu); it != patternGraph.inEnd(pu); ++it) {
                node u = *it;
                if (core1[u] == none && in1[u] == none && out1[u] == none) {
                    if (++in1Neighbors > in2Neighbors) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    /// @ref none is a wildcard on either side.
    bool ruleLabels(node pu, node tv) const {

        if (!nodeLabelled) {
            return true;
        }

        if ((*patternNodeLabels)[pu] == (*targetNodeLabels)[tv] || (*patternNodeLabels)[pu] == none
            || (*targetNodeLabels)[tv] == none) {
            return true;
        }

        return false;
    }

    /// The smallest node in a terminal set, or @ref none if the set is empty. @a positions holds
    /// @ref none for members that have left the set.
    static node smallestMember(const std::vector<node> &members,
                               const std::vector<index> &positions) {
        node smallest = none;
        for (node u : members) {
            // none is the largest representable id, so an empty set falls out of the comparison.
            if (positions[u] != none && u < smallest) {
                smallest = u;
            }
        }
        return smallest;
    }

    /// Removes the members that joined a terminal set after it had size @a mark.
    static void popTail(std::vector<node> &members, std::vector<index> &positions, count mark,
                        count &size) {
        for (index i = mark; i < members.size(); ++i) {
            positions[members[i]] = none;
        }
        size -= members.size() - mark;
        members.resize(mark);
    }

    /**
     * Maps @a pu to @a tv and extends the terminal sets. Returns the record removePair() needs:
     * entries 0 to 3 hold the positions of @a pu and @a tv in the member vectors, or @ref none,
     * and entries 4 to 7 hold the sizes of the member vectors.
     */
    std::array<count, 8> addPair(node pu, node tv) {

        // The sizes go in below, once the two nodes have left their sets but before any neighbour
        // joins one.
        std::array<count, 8> restoreTerminalSets = {in1[pu], out1[pu], in2[tv], out2[tv],
                                                    0,       0,        0,       0};

        core1[pu] = tv;
        core2[tv] = pu;

        if (in1[pu] != none) {
            in1[pu] = none;
            t1in--;
        }
        if (out1[pu] != none) {
            out1[pu] = none;
            t1out--;
        }
        if (in2[tv] != none) {
            in2[tv] = none;
            t2in--;
        }
        if (out2[tv] != none) {
            out2[tv] = none;
            t2out--;
        }

        restoreTerminalSets[4] = membersIn1.size();
        restoreTerminalSets[5] = membersIn2.size();
        restoreTerminalSets[6] = membersOut1.size();
        restoreTerminalSets[7] = membersOut2.size();

        // An undirected graph has in1 = out1 and in2 = out2, so the out-neighbours suffice.
        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it) {
            node u = *it;
            if (core1[u] == none && out1[u] == none) {
                out1[u] = membersOut1.size();
                membersOut1.push_back(u);
                t1out++;

                if (!patternGraph.isDirected()) {
                    in1[u] = membersIn1.size();
                    membersIn1.push_back(u);
                    t1in++;
                }
            }
        }

        for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
            node v = *it;
            if (core2[v] == none && out2[v] == none) {
                out2[v] = membersOut2.size();
                membersOut2.push_back(v);
                t2out++;

                if (!targetGraph.isDirected()) {
                    in2[v] = membersIn2.size();
                    membersIn2.push_back(v);
                    t2in++;
                }
            }
        }

        if (patternGraph.isDirected()) {
            for (auto it = patternGraph.inBegin(pu); it != patternGraph.inEnd(pu); ++it) {
                node u = *it;
                if (core1[u] == none && in1[u] == none) {
                    in1[u] = membersIn1.size();
                    membersIn1.push_back(u);
                    t1in++;
                }
            }
        }

        if (targetGraph.isDirected()) {
            for (auto it = targetGraph.inBegin(tv); it != targetGraph.inEnd(tv); ++it) {
                node v = *it;
                if (core2[v] == none && in2[v] == none) {
                    in2[v] = membersIn2.size();
                    membersIn2.push_back(v);
                    t2in++;
                }
            }
        }

        return restoreTerminalSets;
    }

    /// Undoes addPair() with the record it returned for this pair.
    void removePair(node pu, node tv, const std::array<count, 8> &restoreTerminalSets) {

        core1[pu] = none;
        core2[tv] = none;

        popTail(membersIn1, in1, restoreTerminalSets[4], t1in);
        popTail(membersIn2, in2, restoreTerminalSets[5], t2in);
        popTail(membersOut1, out1, restoreTerminalSets[6], t1out);
        popTail(membersOut2, out2, restoreTerminalSets[7], t2out);

        if (restoreTerminalSets[0] != none) {
            in1[pu] = restoreTerminalSets[0];
            t1in++;
        }
        if (restoreTerminalSets[1] != none) {
            out1[pu] = restoreTerminalSets[1];
            t1out++;
        }
        if (restoreTerminalSets[2] != none) {
            in2[tv] = restoreTerminalSets[2];
            t2in++;
        }
        if (restoreTerminalSets[3] != none) {
            out2[tv] = restoreTerminalSets[3];
            t2out++;
        }
    }

    bool reportMapping() {

        for (index u = 0; u < mapping.size(); ++u) {
            if (patternGraph.hasNode(u)) {
                mapping[u] = core1[u];
            }
        }

        return report(mapping);
    }

    SearchGraph patternGraph;
    SearchGraph targetGraph;

    const std::vector<index> *patternNodeLabels;
    const std::vector<index> *targetNodeLabels;
    bool nodeLabelled;
    bool edgeLabelled;

    SubgraphIsomorphism::Semantics semantics;

    Aux::SignalHandler *handler;

    MatchReporter report;

    /// `core1[u]` is the target node that pattern node u is mapped to, or @ref none.
    std::vector<node> core1;
    /// `core2[v]` is the pattern node mapped to target node v, or @ref none.
    std::vector<node> core2;

    /// Position of each node in the member vectors, or @ref none if it is not in the terminal set.
    std::vector<index> in1, out1, in2, out2;
    /// Current sizes of the four terminal sets.
    count t1in, t1out, t2in, t2out;

    std::vector<node> mapping;

    std::vector<node> membersIn1, membersIn2, membersOut1, membersOut2;
};

} // namespace

VF2::VF2(const Graph &pattern, const Graph &target, Semantics semantics, count maxMatches)
    : SubgraphIsomorphism(pattern, target, semantics, maxMatches) {}

void VF2::run() {
    Aux::SignalHandler handler;
    prepareRun();
    VF2Impl(*pattern, *target, patternNodeLabels, targetNodeLabels, patternEdgeLabels,
            targetEdgeLabels, semantics, handler,
            [this](const Match &match) { return reportMatch(match); })
        .run();
    finishRun();
}

} // namespace NetworKit
