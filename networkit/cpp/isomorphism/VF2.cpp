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
    /**
     * @param pattern Snapshot of the pattern, built with the adjacency matrix.
     * @param target Snapshot of the target, built without it.
     * @param patternNodeLabels Empty when the search is unlabelled.
     * @param targetNodeLabels Empty when the search is unlabelled.
     * @param semantics Whether matches must be induced.
     * @param handler Polled so a long search can be stopped with CTRL+C.
     * @param report Where complete mappings are reported.
     */
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

    /**
     * Search for every match and report each one. Initialize core1/core2 and in1/out1/in2/out2.
     * Handle trivial cases immediately and call match(0) to start the recursion.
     */
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
    /**
     * One level of the depth-first search: extend a mapping of @a depth pairs by one more.
     *
     * @param depth Current search depth.
     * @return false if the whole search must stop, true otherwise.
     */
    bool match(count depth) {

        // If depth equals the number of pattern nodes, all pattern nodes are mapped
        if (depth == patternGraph.numberOfNodes()) {
            return reportMapping();
        }

        index cursor = 0;
        node pu = none;
        node tv = none;
        bool continueSearch;

        // Iterate over all candidate pairs and if candidate pair is feasible, add pair and call
        // match(depth + 1)
        while (nextCandidatePair(cursor, pu, tv)) {
            handler->assureRunning();
            if (feasible(pu, tv)) {
                const std::array<count, 8> restoreTerminalSets = addPair(pu, tv);
                continueSearch = match(depth + 1);
                // Remove pair independent of outcome and abort search if it must be stopped
                removePair(pu, tv, restoreTerminalSets);
                if (!continueSearch) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Produce the next candidate pair to try at this depth.
     *
     * @param depth Current search depth.
     * @param cursor In/out: where the previous call stopped, so iteration can resume.
     * @param pu Out: the pattern node to map.
     * @param tv Out: the target node to try for it.
     * @return false when the candidates at this depth are exhausted.
     */
    bool nextCandidatePair(index &cursor, node &pu, node &tv) const {

        if (t1out != 0 && t2out != 0) {

            // Smallest pattern node still in out1. t1out is nonzero, so one exists.
            pu = smallestMember(membersOut1, out1);
            if (pu == none) {
                return false;
            }
            // Pair it with every target node still in out2
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

            // Smallest pattern node still in in1. t1in is nonzero, so one exists.
            pu = smallestMember(membersIn1, in1);
            if (pu == none) {
                return false;
            }
            // Pair it with every target node still in in2
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

            // Find smallest unmapped pattern node
            for (node u = 0; u < core1.size(); ++u) {
                if (patternGraph.hasNode(u) && core1[u] == none) {
                    pu = u;
                    break;
                }
            }
            // Pair it with every unmapped target node
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

    /**
     * @param pu The pattern node.
     * @param tv The target node.
     * @return true if the pair @a pu, @a tv may be added to the mapping.
     */

    bool feasible(node pu, node tv) const {

        return ruleSuccessors(pu, tv) && rulePredecessors(pu, tv) && ruleTerminalCounts(pu, tv)
               && ruleNewCounts(pu, tv) && ruleLabels(pu, tv);
    }

    /**
     * Consistency check for out-edges. For every out-neighbour of @a pu that is already mapped, the
     * target must contain the corresponding edge out of @a tv.
     *
     * @param pu The pattern node.
     * @param tv The target node.
     * @return true if the out-edges of @a pu, @a tv are consistent.
     */
    bool ruleSuccessors(node pu, node tv) const {

        // Check for every mapped out-neighbor of pu if the target has the corresponding edge out of
        // tv
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

        // Under Semantics::INDUCED: Check for every mapped out-neighbor of tv if the pattern has
        // the corresponding edge out of pu
        if (semantics == SubgraphIsomorphism::Semantics::INDUCED) {
            for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
                node v = *it;
                if (core2[v] != none) {
                    if (!patternGraph.hasEdge(pu, core2[v])) {
                        return false;
                    }
                    if (edgeLabelled) {
                        if (!(patternGraph.edgeLabel(tv, v) == targetGraph.edgeLabel(pu, core2[v])
                              || patternGraph.edgeLabel(tv, v) == none
                              || targetGraph.edgeLabel(pu, core2[v]) == none)) {
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    /**
     * Consistency check for in-edges. For every in-neighbour of @a pu that is already mapped, the
     * target must contain the corresponding edge into @a tv. The mirror image of @ref
     * ruleSuccessors().
     *
     * @param pu The pattern node.
     * @param tv The target node.
     * @return true if the in-edges of @a pu, @a tv are consistent.
     */
    bool rulePredecessors(node pu, node tv) const {

        // If undirected, in-neighbors=out-neighbors, return true immediately
        if (!patternGraph.isDirected()) {
            return true;
        }

        // Check for every mapped in-neighbor of pu if the target has the corresponding edge into tv
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

        // Under Semantics::INDUCED: Check for every mapped in-neighbor of tv if the pattern has the
        // corresponding edge into pu
        if (semantics == SubgraphIsomorphism::Semantics::INDUCED) {
            for (auto it = targetGraph.inBegin(tv); it != targetGraph.inEnd(tv); ++it) {
                node v = *it;
                if (core2[v] != none) {
                    if (!patternGraph.hasEdge(core2[v], pu)) {
                        return false;
                    }
                    if (edgeLabelled) {
                        if (!(patternGraph.edgeLabel(v, tv) == targetGraph.edgeLabel(core2[v], pu)
                              || patternGraph.edgeLabel(v, tv) == none
                              || targetGraph.edgeLabel(core2[v], pu) == none)) {
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    /**
     * One-step look-ahead on the terminal sets. Count unmapped, out- and in-terminal neighbors of
     * @a pu and @a tv. Return false if the pattern count exceeds the target count for either.
     *
     * @param pu The pattern node.
     * @param tv The target node.
     * @return true if one-step look-ahead on the terminal sets of @a pu, @a tv passes.
     */
    bool ruleTerminalCounts(node pu, node tv) const {

        count in1Neighbors = 0;
        count out1Neighbors = 0;
        count in2Neighbors = 0;
        count out2Neighbors = 0;

        // Count unmapped, out-terminal neighbors of tv and pu
        for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
            node v = *it;
            if (core2[v] == none && out2[v] != none) {
                out2Neighbors++;
            }
        }
        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it) {
            node u = *it;
            if (core1[u] == none && out1[u] != none) {
                // Return false if pu has more such neighbors than tv
                if (++out1Neighbors > out2Neighbors) {
                    return false;
                }
            }
        }

        // If directed, do the same for unmapped, in-terminal neighbors
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

    /**
     * Two-step look-ahead on the terminal sets. Count unmapped, out- and in-neighbors of @a pu and
     * @a tv that are not part of any terminal set. Return false if the pattern count exceeds the
     * target count for either.
     *
     * @param pu The pattern node.
     * @param tv The target node.
     * @return true if two-step look-ahead on the terminal sets of @a pu, @a tv passes.
     */
    bool ruleNewCounts(node pu, node tv) const {

        // Semantics::MONOMORPHISM allows extra target edges, return true immediately
        if (semantics == SubgraphIsomorphism::Semantics::MONOMORPHISM) {
            return true;
        }

        count in1Neighbors = 0;
        count out1Neighbors = 0;
        count in2Neighbors = 0;
        count out2Neighbors = 0;

        // Count unmapped, non-terminal out-neighbors of tv and pu
        for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
            node v = *it;
            if (core2[v] == none && in2[v] == none && out2[v] == none) {
                out2Neighbors++;
            }
        }
        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it) {
            node u = *it;
            if (core1[u] == none && in1[u] == none && out1[u] == none) {
                // Return false if pu has more such neighbors than tv
                if (++out1Neighbors > out2Neighbors) {
                    return false;
                }
            }
        }

        // If directed, do the same for in-neighbors
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

    /**
     * Consistency check for labels. Return true immediately if the search is unlabelled. Otherwise
     * the labels of @a pu and @a tv must be equal, except that @ref none on either side is a
     * wildcard that matches anything.
     *
     * @param pu The pattern node.
     * @param tv The target node.
     * @return true if the labels of @a pu, @a tv are consistent.
     */
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

    /**
     * The smallest node still in a terminal set, or @ref none when the set holds nothing.
     *
     * @param members The member vector of the terminal set.
     * @param positions The corresponding in1/out1/in2/out2 vector storing index of each node in @a
     * members. An entry of @ref none means that the node has been mapped and is no longer part of
     * the terminal set.
     * @return the smallest node still in a terminal set, or @ref none when the set holds nothing.
     */
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

    /**
     * Drop everything a terminal set gained since it had size @a mark and set the terminal set
     * member vector positions of the removed nodes to @ref none.
     *
     * @param members The member vector of the terminal set.
     * @param positions The corresponding in1/out1/in2/out2 vector storing index of each node in @a
     * members.
     * @param mark The index from which to remove the nodes.
     */
    static void popTail(std::vector<node> &members, std::vector<index> &positions, count mark,
                        count &size) {
        for (index i = mark; i < members.size(); ++i) {
            positions[members[i]] = none;
        }
        size -= members.size() - mark;
        members.resize(mark);
    }

    /**
     * Add (@a pu, @a tv) to the mapping and update the four terminal sets.
     *
     * @param pu The pattern node.
     * @param tv The target node.
     * @return The record @ref removePair() needs to undo this call. Entries 0 to 3 hold the
     * positions @a pu and @a tv had in the four member vectors, or @ref none where the node was
     * not in that set. Entries 4 to 7 hold the sizes the four member vectors had beforehand.
     */
    std::array<count, 8> addPair(node pu, node tv) {

        // Where pu and tv sit in the member vectors right now. The sizes go in below, once the
        // two nodes have left their sets but before any neighbour joins one.
        std::array<count, 8> restoreTerminalSets = {in1[pu], out1[pu], in2[tv], out2[tv],
                                                    0,       0,        0,       0};

        // Map pu and tv onto each other
        core1[pu] = tv;
        core2[tv] = pu;

        // Reset the positions that pu and tv have in the terminal set member vectors to none
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

        // Store the current size of the terminal set member vectors
        restoreTerminalSets[4] = membersIn1.size();
        restoreTerminalSets[5] = membersIn2.size();
        restoreTerminalSets[6] = membersOut1.size();
        restoreTerminalSets[7] = membersOut2.size();

        // Unmapped neighbors of pu and tv are added to the respective terminal sets
        // If undirected, iterating over out-neighbors is sufficient, because in1=out1 and in2=out2
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

        // If directed, iterate over inNeighbors separately
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

    /**
     * Undo @ref addPair() exactly.
     *
     * @param pu The pattern node.
     * @param tv The target node.
     * @param restoreTerminalSets The record @ref addPair() returned for this very pair.
     */
    void removePair(node pu, node tv, const std::array<count, 8> &restoreTerminalSets) {

        // Unmap pu and tv
        core1[pu] = none;
        core2[tv] = none;

        // Drop everything the four terminal sets gained in addPair(pu, tv) and set the terminal set
        // member vector positions of removed nodes to none
        popTail(membersIn1, in1, restoreTerminalSets[4], t1in);
        popTail(membersIn2, in2, restoreTerminalSets[5], t2in);
        popTail(membersOut1, out1, restoreTerminalSets[6], t1out);
        popTail(membersOut2, out2, restoreTerminalSets[7], t2out);

        // If pu or tv were part of any terminal sets before addPair(pu, tv), add them back in and
        // restore the positions they had in the terminal set member vectors before addPair(pu, tv)
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

    /**
     * Hand a complete mapping over. Copy core1 into 'mapping' for the pattern nodes that exist,
     * then report the mapping.
     */
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

    /// Signal handler used to abort the search on interruption.
    Aux::SignalHandler *handler;

    MatchReporter report;

    /// core1[patternNode] = target node it is mapped to, or `none`.
    std::vector<node> core1;
    /// core2[targetNode] = pattern node mapped onto it, or `none`.
    std::vector<node> core2;

    /// Index at which each node can be found in the terminal set member vectors; none means "not in
    /// it".
    std::vector<index> in1, out1, in2, out2;
    /// Current sizes of the four terminal sets.
    count t1in, t1out, t2in, t2out;

    /// Reused buffer handed to the reporter, so a match costs no allocation.
    std::vector<node> mapping;

    /// Member vectors for in and out terminal sets.
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
