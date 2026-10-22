#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/isomorphism/VF2.hpp>

#include "MatchReporter.hpp"
#include "SearchGraph.hpp"

namespace NetworKit {

namespace {

using IsomorphismDetails::labelsCompatible;
using IsomorphismDetails::MatchReporter;
using IsomorphismDetails::nodeLabelsCompatible;
using IsomorphismDetails::SearchGraph;

/// The unmapped nodes of one graph that are adjacent to mapped nodes in one direction.
class TerminalSet {

public:
    /// What restore() needs to undo leave() and the extend() calls after it.
    struct Checkpoint {
        bool wasMember;
        count mark;
    };

    explicit TerminalSet(count upperNodeIdBound) : isMember(upperNodeIdBound, false) {}

    bool empty() const noexcept { return size == 0; }

    bool contains(node u) const noexcept {
        // An undirected search gives its in-sets no storage, so it must never query them.
        assert(u < isMember.size());
        return isMember[u];
    }

    /// The nodes in the order they joined. Nodes that left the set may stay listed, so callers
    /// skip the entries that contains() rejects.
    const std::vector<node> &joined() const noexcept { return members; }

    /// The smallest node in the set, or @ref none if the set is empty.
    node smallest() const {
        node result = none;
        for (const node u : members) {
            // none is the largest representable id, so an empty set falls out of the comparison.
            if (contains(u) && u < result)
                result = u;
        }
        return result;
    }

    /// Removes @a u, which has just been mapped, from the set.
    Checkpoint leave(node u) {
        const Checkpoint checkpoint{contains(u), members.size()};
        if (checkpoint.wasMember) {
            isMember[u] = false;
            --size;
        }
        return checkpoint;
    }

    /// Adds the nodes in [begin, end) that are unmapped in @a core and not yet in the set.
    void extend(const node *begin, const node *end, const std::vector<node> &core) {
        for (const node *it = begin; it != end; ++it) {
            const node u = *it;
            if (core[u] == none && !isMember[u]) {
                isMember[u] = true;
                members.push_back(u);
                ++size;
            }
        }
    }

    /// Removes the nodes that joined after leave(@a u) returned @a checkpoint, then puts @a u back.
    void restore(node u, Checkpoint checkpoint) {
        for (index i = checkpoint.mark; i < members.size(); ++i)
            isMember[members[i]] = false;
        size -= members.size() - checkpoint.mark;
        members.resize(checkpoint.mark);

        if (checkpoint.wasMember) {
            isMember[u] = true;
            ++size;
        }
    }

private:
    /// Whether each node is in the set. One byte per node keeps the target-side sets small.
    std::vector<uint8_t> isMember;
    std::vector<node> members;
    count size = 0;
};

class VF2Impl {

public:
    VF2Impl(const Graph &pattern, const Graph &target, const std::vector<index> &patternNodeLabels,
            const std::vector<index> &targetNodeLabels, const std::vector<index> &patternEdgeLabels,
            const std::vector<index> &targetEdgeLabels, SubgraphIsomorphism::Semantics semantics,
            Aux::SignalHandler &handler, MatchReporter report)
        : patternGraph(pattern, /* buildMatrix = */ true, patternEdgeLabels),
          targetGraph(target, /* buildMatrix = */ false, targetEdgeLabels),
          patternNodeLabels(&patternNodeLabels), targetNodeLabels(&targetNodeLabels),
          induced(semantics == SubgraphIsomorphism::Semantics::INDUCED), handler(&handler),
          report(std::move(report)), core1(patternGraph.upperNodeIdBound(), none),
          core2(targetGraph.upperNodeIdBound(), none),
          // An undirected graph uses only the out-sets, so its in-sets get no storage.
          in1(directed() ? patternGraph.upperNodeIdBound() : 0),
          out1(patternGraph.upperNodeIdBound()),
          in2(directed() ? targetGraph.upperNodeIdBound() : 0),
          out2(targetGraph.upperNodeIdBound()) {
        IsomorphismDetails::refuseLostLabels(patternGraph, targetGraph, "VF2");
    }

    void run() {
        if (patternGraph.numberOfNodes() > targetGraph.numberOfNodes()
            || patternGraph.maxInDegree() > targetGraph.maxInDegree()
            || patternGraph.maxOutDegree() > targetGraph.maxOutDegree())
            return;

        search();
    }

private:
    /// One level of the depth-first search: a pattern node and the state of its candidates.
    struct Frame {
        node pu;
        /// The terminal set that holds the candidate images of `pu`, or nullptr if every unmapped
        /// target node is a candidate.
        const TerminalSet *candidates;
        /// The position of the next candidate in its list.
        index cursor;
        TerminalSet::Checkpoint in1Saved, out1Saved, in2Saved, out2Saved;
    };

    /// Both graphs agree on directedness, see SubgraphIsomorphism::validateInput().
    bool directed() const noexcept { return patternGraph.isDirected(); }

    /// The depth-first search. Each frame maps one pattern node, and the frames live on the heap,
    /// so a large pattern cannot overflow the call stack.
    void search() {
        const count depthOfMatch = patternGraph.numberOfNodes();

        // A pattern without nodes has one match, the empty mapping.
        if (depthOfMatch == 0) {
            report(core1);
            return;
        }

        std::vector<Frame> stack;
        stack.reserve(depthOfMatch);
        pushFrame(stack);

        // assureRunning() is an out-of-line call, and polling it for every candidate costs 4-7%.
        constexpr count pollInterval = 1024;
        count candidatesTried = 0;

        while (!stack.empty()) {
            Frame &frame = stack.back();
            // A mapped pattern node means that the search has come back from a deeper frame.
            if (core1[frame.pu] != none)
                unmapPair(frame);

            const node tv = nextCandidate(frame);
            if (tv == none) {
                stack.pop_back();
                continue;
            }

            if (++candidatesTried % pollInterval == 0)
                handler->assureRunning();
            if (!feasible(frame.pu, tv))
                continue;

            if (stack.size() < depthOfMatch) {
                mapPair(frame, tv);
                pushFrame(stack);
                continue;
            }

            // The last pair completes a match. No frame reads the terminal sets after it, so
            // only core1 needs the pair. It then holds the image of every pattern node.
            core1[frame.pu] = tv;
            const bool continueSearch = report(core1);
            core1[frame.pu] = none;
            if (!continueSearch)
                return;
        }
    }

    /// Pushes the frame of the next pattern node, unless the terminal sets show a dead end.
    void pushFrame(std::vector<Frame> &stack) const {
        // A pattern node in a terminal set needs an image in the matching target set.
        if ((!out1.empty() && out2.empty()) || (!in1.empty() && in2.empty()))
            return;

        Frame frame{};
        if (!out1.empty()) {
            frame.pu = out1.smallest();
            frame.candidates = &out2;
        } else if (!in1.empty()) {
            frame.pu = in1.smallest();
            frame.candidates = &in2;
        } else {
            // The frame keeps no candidate set, so every unmapped target node is a candidate.
            while (!patternGraph.hasNode(frame.pu) || core1[frame.pu] != none)
                ++frame.pu;
        }
        stack.push_back(frame);
    }

    /// Returns the next candidate image of `frame.pu` and advances the cursor past it. Returns
    /// @ref none once the candidates are exhausted.
    node nextCandidate(Frame &frame) const {
        if (frame.candidates == nullptr) {
            while (frame.cursor < core2.size()) {
                const node tv = frame.cursor++;
                if (targetGraph.hasNode(tv) && core2[tv] == none)
                    return tv;
            }
            return none;
        }

        // mapPair() appends to the member list, and unmapPair() truncates it again. Every pair
        // mapped at this depth or below is unmapped before control returns here, so the list is the
        // same at every visit of this frame.
        const std::vector<node> &members = frame.candidates->joined();
        while (frame.cursor < members.size()) {
            const node tv = members[frame.cursor++];
            if (frame.candidates->contains(tv))
                return tv;
        }
        return none;
    }

    /// Maps `frame.pu` to @a tv and extends the terminal sets. unmapPair() undoes both.
    void mapPair(Frame &frame, node tv) {
        const node pu = frame.pu;
        core1[pu] = tv;
        core2[tv] = pu;

        frame.out1Saved = out1.leave(pu);
        frame.out2Saved = out2.leave(tv);
        out1.extend(patternGraph.outBegin(pu), patternGraph.outEnd(pu), core1);
        out2.extend(targetGraph.outBegin(tv), targetGraph.outEnd(tv), core2);

        if (directed()) {
            frame.in1Saved = in1.leave(pu);
            frame.in2Saved = in2.leave(tv);
            in1.extend(patternGraph.inBegin(pu), patternGraph.inEnd(pu), core1);
            in2.extend(targetGraph.inBegin(tv), targetGraph.inEnd(tv), core2);
        }
    }

    void unmapPair(const Frame &frame) {
        const node pu = frame.pu;
        const node tv = core1[pu];
        core1[pu] = none;
        core2[tv] = none;

        out1.restore(pu, frame.out1Saved);
        out2.restore(tv, frame.out2Saved);
        if (directed()) {
            in1.restore(pu, frame.in1Saved);
            in2.restore(tv, frame.in2Saved);
        }
    }

    bool feasible(node pu, node tv) const {
        // The node labels take constant time to check, so they go first.
        return nodeLabelsCompatible(*patternNodeLabels, *targetNodeLabels, pu, tv)
               && ruleSuccessors(pu, tv) && rulePredecessors(pu, tv)
               && (!patternGraph.hasEdgeLabels() || ruleEdgeLabels(pu, tv))
               && ruleTerminalCounts(pu, tv) && ruleNewCounts(pu, tv);
    }

    /// Every mapped out-neighbour of @a pu needs a corresponding out-edge of @a tv. Under INDUCED
    /// semantics, the same holds in reverse.
    bool ruleSuccessors(node pu, node tv) const {
        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it) {
            const node u = *it;
            if (core1[u] != none && !targetGraph.hasEdge(tv, core1[u]))
                return false;
        }

        // ruleEdgeLabels() compares the labels of the edges that both graphs contain.
        if (induced) {
            for (auto it = targetGraph.outBegin(tv); it != targetGraph.outEnd(tv); ++it) {
                const node v = *it;
                if (core2[v] != none && !patternGraph.hasEdge(pu, core2[v]))
                    return false;
            }
        }

        return true;
    }

    /// The mirror image of ruleSuccessors() for in-edges.
    bool rulePredecessors(node pu, node tv) const {
        if (!directed())
            return true;

        for (auto it = patternGraph.inBegin(pu); it != patternGraph.inEnd(pu); ++it) {
            const node u = *it;
            if (core1[u] != none && !targetGraph.hasEdge(core1[u], tv))
                return false;
        }

        if (induced) {
            for (auto it = targetGraph.inBegin(tv); it != targetGraph.inEnd(tv); ++it) {
                const node v = *it;
                if (core2[v] != none && !patternGraph.hasEdge(core2[v], pu))
                    return false;
            }
        }

        return true;
    }

    /// The images of the edges between @a pu and the mapped pattern nodes need compatible labels.
    /// ruleSuccessors() and rulePredecessors() have checked that the images exist. A separate pass
    /// keeps the label loads out of their loops when the graphs carry no edge labels.
    bool ruleEdgeLabels(node pu, node tv) const {
        const index *label = patternGraph.outLabelBegin(pu);
        for (auto it = patternGraph.outBegin(pu); it != patternGraph.outEnd(pu); ++it, ++label) {
            const node u = *it;
            if (core1[u] != none && !labelsCompatible(*label, targetGraph.edgeLabel(tv, core1[u])))
                return false;
        }

        if (!directed())
            return true;

        label = patternGraph.inLabelBegin(pu);
        for (auto it = patternGraph.inBegin(pu); it != patternGraph.inEnd(pu); ++it, ++label) {
            const node u = *it;
            if (core1[u] != none && !labelsCompatible(*label, targetGraph.edgeLabel(core1[u], tv)))
                return false;
        }

        return true;
    }

    /// One-step look-ahead: @a pu must not have more unmapped neighbours in the terminal sets than
    /// @a tv has.
    bool ruleTerminalCounts(node pu, node tv) const {
        // A terminal set holds only unmapped nodes, so membership alone decides.
        const auto inOut1 = [this](node u) { return out1.contains(u); };
        const auto inOut2 = [this](node v) { return out2.contains(v); };
        if (!enoughInTarget(patternGraph.outBegin(pu), patternGraph.outEnd(pu), inOut1,
                            targetGraph.outBegin(tv), targetGraph.outEnd(tv), inOut2))
            return false;

        if (!directed())
            return true;

        const auto inIn1 = [this](node u) { return in1.contains(u); };
        const auto inIn2 = [this](node v) { return in2.contains(v); };
        return enoughInTarget(patternGraph.inBegin(pu), patternGraph.inEnd(pu), inIn1,
                              targetGraph.inBegin(tv), targetGraph.inEnd(tv), inIn2);
    }

    /// Two-step look-ahead: the same count for unmapped neighbours outside all terminal sets.
    bool ruleNewCounts(node pu, node tv) const {
        // MONOMORPHISM allows extra target edges, so this rule does not apply.
        if (!induced)
            return true;

        // An undirected graph has no in-sets to look up.
        const auto new1 = [this](node u) {
            return core1[u] == none && !out1.contains(u) && !(directed() && in1.contains(u));
        };
        const auto new2 = [this](node v) {
            return core2[v] == none && !out2.contains(v) && !(directed() && in2.contains(v));
        };
        return enoughInTarget(patternGraph.outBegin(pu), patternGraph.outEnd(pu), new1,
                              targetGraph.outBegin(tv), targetGraph.outEnd(tv), new2)
               && (!directed()
                   || enoughInTarget(patternGraph.inBegin(pu), patternGraph.inEnd(pu), new1,
                                     targetGraph.inBegin(tv), targetGraph.inEnd(tv), new2));
    }

    /// Whether [targetBegin, targetEnd) holds at least as many nodes that pass @a targetTest as
    /// [patternBegin, patternEnd) holds nodes that pass @a patternTest.
    template <typename PatternTest, typename TargetTest>
    static bool enoughInTarget(const node *patternBegin, const node *patternEnd,
                               PatternTest patternTest, const node *targetBegin,
                               const node *targetEnd, TargetTest targetTest) {
        // The pattern side is small, and its count lets the target scan stop early.
        auto needed = std::count_if(patternBegin, patternEnd, patternTest);
        for (const node *it = targetBegin; needed != 0 && it != targetEnd; ++it) {
            if (targetTest(*it))
                --needed;
        }
        return needed == 0;
    }

    SearchGraph patternGraph;
    SearchGraph targetGraph;

    const std::vector<index> *patternNodeLabels;
    const std::vector<index> *targetNodeLabels;

    /// True for SubgraphIsomorphism::Semantics::INDUCED, false for MONOMORPHISM.
    bool induced;

    Aux::SignalHandler *handler;

    MatchReporter report;

    /// `core1[u]` is the target node that pattern node u is mapped to, or @ref none.
    std::vector<node> core1;
    /// `core2[v]` is the pattern node mapped to target node v, or @ref none.
    std::vector<node> core2;

    TerminalSet in1, out1, in2, out2;
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
