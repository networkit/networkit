#include <stdexcept>
#include <vector>

#include <tlx/unused.hpp>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/isomorphism/VF3.hpp>

#include "MatchReporter.hpp"
#include "SearchGraph.hpp"

namespace NetworKit {

namespace {

using IsomorphismDetails::MatchReporter;
using IsomorphismDetails::SearchGraph;

/**
 * The VF3 search. Not implemented yet; VF3::run() throws before constructing it.
 *
 * The mapping is stored in both directions, `core1` and `core2`, as in VF2. Everything else is
 * precomputed before the search and does not change afterwards:
 *
 * - `patternClass` and `targetClass` assign every node a class derived from its label, and
 *   `targetClassSize` counts the target nodes of each class.
 * - `order` is the fixed matching order, and `orderParent[i]` is the position of an earlier
 *   neighbour of `order[i]`, or `none` if position @a i starts a new component. The candidates
 *   for position @a i are drawn from the neighbours of the parent's image.
 * - `feasibilitySets[depth][cls]` is the number of nodes of class `cls` that the pattern still
 *   needs from `depth` on.
 */
class VF3Impl {

public:
    /**
     * @param pattern The pattern graph.
     * @param target The target graph.
     * @param patternNodeLabels Empty if the search is unlabelled; then there is one class.
     * @param targetNodeLabels Empty if the search is unlabelled.
     * @param semantics Whether matches must be induced.
     * @param handler Signal handler used to abort the search on interruption.
     * @param report Receives every complete mapping.
     */
    VF3Impl(const Graph &pattern, const Graph &target, const std::vector<index> &patternNodeLabels,
            const std::vector<index> &targetNodeLabels, SubgraphIsomorphism::Semantics semantics,
            Aux::SignalHandler &handler, MatchReporter report)
        : patternGraph(pattern, /* buildMatrix = */ true),
          targetGraph(target, /* buildMatrix = */ false), patternNodeLabels(&patternNodeLabels),
          targetNodeLabels(&targetNodeLabels), nodeLabelled(!patternNodeLabels.empty()),
          semantics(semantics), handler(&handler), report(std::move(report)), numberOfClasses(0) {}

    /**
     * Searches for every match and reports each one.
     *
     * TODO: implement.
     *  1. classifyNodes(), then computeNodeOrder(), then precomputeFeasibilitySets(). Each step
     *     needs the previous one.
     *  2. Size core1, core2 and the reusable `mapping` buffer, filling the cores with `none`.
     *  3. Return early if any class needs more pattern nodes than the target has.
     *  4. Call match(0).
     */
    void run() {
        // TODO: remove once implemented.
        tlx::unused(patternGraph, targetGraph, patternNodeLabels, targetNodeLabels, nodeLabelled,
                    semantics, handler, report, patternClass, targetClass, targetClassSize,
                    numberOfClasses, order, orderParent, feasibilitySets, core1, core2, mapping);
        throw std::logic_error("VF3Impl::run() is not implemented yet");
    }

private:
    /**
     * Groups the nodes of both graphs into classes.
     *
     * TODO: implement. Map each distinct label to a small dense class id and fill
     * patternClass/targetClass with it, counting the target's classes into targetClassSize. If the
     * search is unlabelled, every node goes into class 0.
     */
    void classifyNodes() {
        throw std::logic_error("VF3Impl::classifyNodes() is not implemented yet");
    }

    /**
     * Computes the fixed order in which pattern nodes are mapped.
     *
     * TODO: implement. Greedily append the pattern node that scores best on, in order of
     * importance: most edges into the nodes already ordered, then lowest classProbability(), then
     * highest degree. Record in orderParent[i] the position of an already ordered neighbour, or
     * `none` if there is none. Section 4 of the paper.
     */
    void computeNodeOrder() {
        throw std::logic_error("VF3Impl::computeNodeOrder() is not implemented yet");
    }

    /**
     * Returns the probability that a random target node belongs to class @a cls.
     *
     * TODO: implement as targetClassSize[cls] divided by the number of target nodes.
     */
    double classProbability(index cls) const {
        tlx::unused(cls);
        throw std::logic_error("VF3Impl::classProbability() is not implemented yet");
    }

    /**
     * Precomputes, per depth and per class, how many nodes the rest of the pattern still needs.
     *
     * TODO: implement. Walk `order` backwards accumulating the class counts, so that
     * feasibilitySets[depth][cls] is the number of nodes of class cls that the pattern requires
     * from depth on. ruleClassCounts() compares this with what the target still has available.
     */
    void precomputeFeasibilitySets() {
        throw std::logic_error("VF3Impl::precomputeFeasibilitySets() is not implemented yet");
    }

    /**
     * One level of the depth-first search: maps the pattern node at position @a depth in `order`.
     *
     * TODO: implement.
     *  1. If @a depth equals the number of pattern nodes, report the mapping and return the
     *     reporter's answer.
     *  2. Otherwise, for each candidate of this depth that passes feasible(), call addPair(),
     *     recurse and call removePair(). Return false as soon as a recursion returns false.
     *
     * @param depth Current position in `order`.
     * @return false if the whole search must stop, true otherwise.
     */
    bool match(count depth) {
        tlx::unused(depth);
        throw std::logic_error("VF3Impl::match() is not implemented yet");
    }

    /**
     * Collects the target nodes to try at position @a depth.
     *
     * TODO: implement. If orderParent[depth] is set, the candidates are the neighbours of the
     * parent's image. If it is `none`, every unmapped target node of the right class is a
     * candidate.
     *
     * @param depth Current position in `order`.
     * @param out Filled with the candidates; cleared first.
     */
    void candidatesFor(count depth, std::vector<node> &out) const {
        tlx::unused(depth, out);
        throw std::logic_error("VF3Impl::candidatesFor() is not implemented yet");
    }

    /**
     * Returns whether the pair (@a pu, @a tv) may be added to the mapping at @a depth.
     *
     * TODO: implement by calling ruleLabels(), ruleEdges() and ruleClassCounts(), cheapest first,
     * and returning false at the first failure.
     */
    bool feasible(node pu, node tv, count depth) const {
        tlx::unused(pu, tv, depth);
        throw std::logic_error("VF3Impl::feasible() is not implemented yet");
    }

    /**
     * Consistency check against the mapped nodes.
     *
     * TODO: implement. Every edge between @a pu and a mapped pattern node must have a counterpart
     * between @a tv and the image of that node, in both directions for a directed graph. Under
     * Semantics::INDUCED, a target edge between @a tv and a mapped node without a corresponding
     * pattern edge also rejects the pair.
     */
    bool ruleEdges(node pu, node tv, count depth) const {
        tlx::unused(pu, tv, depth);
        throw std::logic_error("VF3Impl::ruleEdges() is not implemented yet");
    }

    /**
     * Look-ahead on the class counts.
     *
     * TODO: implement. For each class, compare feasibilitySets[depth] with the number of unmapped
     * target nodes of that class that are still reachable. Reject the pair if the pattern needs
     * more nodes of some class than remain.
     */
    bool ruleClassCounts(node pu, node tv, count depth) const {
        tlx::unused(pu, tv, depth);
        throw std::logic_error("VF3Impl::ruleClassCounts() is not implemented yet");
    }

    /**
     * Consistency check for labels.
     *
     * TODO: implement. Return true immediately if the search is unlabelled. Otherwise, the
     * classes must agree, with @ref none acting as a wildcard on either side.
     */
    bool ruleLabels(node pu, node tv) const {
        tlx::unused(pu, tv);
        throw std::logic_error("VF3Impl::ruleLabels() is not implemented yet");
    }

    /**
     * Adds (@a pu, @a tv) to the mapping.
     *
     * TODO: implement. Set core1[pu] = tv and core2[tv] = pu, and update the per-class
     * availability counters that ruleClassCounts() reads.
     */
    void addPair(node pu, node tv, count depth) {
        tlx::unused(pu, tv, depth);
        throw std::logic_error("VF3Impl::addPair() is not implemented yet");
    }

    /**
     * Undoes @ref addPair().
     *
     * TODO: implement. Reset core1[pu] and core2[tv] to `none` and restore the counters that
     * addPair() changed.
     */
    void removePair(node pu, node tv, count depth) {
        tlx::unused(pu, tv, depth);
        throw std::logic_error("VF3Impl::removePair() is not implemented yet");
    }

    /**
     * Reports a complete mapping.
     *
     * TODO: implement. Copy core1 into `mapping`, which is indexed by pattern node, and return
     * report(mapping).
     */
    bool reportMapping() { throw std::logic_error("VF3Impl::reportMapping() is not implemented"); }

    SearchGraph patternGraph;
    SearchGraph targetGraph;

    const std::vector<index> *patternNodeLabels;
    const std::vector<index> *targetNodeLabels;
    bool nodeLabelled;

    SubgraphIsomorphism::Semantics semantics;

    /// Signal handler used to abort the search on interruption.
    Aux::SignalHandler *handler;

    MatchReporter report;

    /// Class of every node, derived from its label. All zero if the search is unlabelled.
    std::vector<index> patternClass, targetClass;
    /// Number of target nodes in each class.
    std::vector<count> targetClassSize;
    count numberOfClasses;

    /// Fixed matching order: order[i] is the pattern node mapped at depth i.
    std::vector<node> order;
    /// orderParent[i] is the position in `order` that position i attaches to, or `none`.
    std::vector<index> orderParent;

    /// feasibilitySets[depth][cls] = nodes of class cls the pattern still needs from `depth` on.
    std::vector<std::vector<count>> feasibilitySets;

    /// core1[patternNode] = target node it is mapped to, or `none`.
    std::vector<node> core1;
    /// core2[targetNode] = pattern node mapped onto it, or `none`.
    std::vector<node> core2;

    /// Reused buffer for reported matches.
    std::vector<node> mapping;
};

} // namespace

VF3::VF3(const Graph &pattern, const Graph &target, Semantics semantics, count maxMatches)
    : SubgraphIsomorphism(pattern, target, semantics, maxMatches) {}

void VF3::run() {
    // Refuse unsupported input before reporting the missing search; the tests tell the two apart.
    if (isEdgeLabelled())
        throw std::runtime_error("VF3 does not support edge labels - see "
                                 "SubgraphIsomorphism::setEdgeLabels()");

    prepareRun();

    // TODO: once VF3Impl is implemented, replace the throw below by
    //
    //     Aux::SignalHandler handler;
    //     VF3Impl(*pattern, *target, patternNodeLabels, targetNodeLabels, semantics, handler,
    //             [this](const Match &match) { return reportMatch(match); }).run();
    //     finishRun();
    throw std::logic_error("VF3 is not implemented yet - see VF3Impl in VF3.cpp");
}

} // namespace NetworKit
