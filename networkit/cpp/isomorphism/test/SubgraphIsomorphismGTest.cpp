/*
 * SubgraphIsomorphismGTest.cpp
 *
 * Tests the base class itself - the run() protocol, the two callback forms, the worker
 * accumulation a parallel algorithm is expected to do, and the interrupt policy.
 *
 * None of this is about any particular search algorithm. It is about the contract the four of
 * them share, which is why it is exercised here through two stand-ins rather than through VF2 or
 * RI: ReferenceSubgraphIsomorphism for the sequential path and MultiWorkerReporter for the
 * parallel one. Both report a precomputed match set, so a failure here is always the base class.
 */

#include <atomic>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <vector>

#include <omp.h>
#include <gtest/gtest.h>

#include <networkit/GlobalState.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/ParallelRI.hpp>
#include <networkit/isomorphism/RI.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>
#include <networkit/isomorphism/VF2.hpp>
#include <networkit/isomorphism/VF3.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"

namespace NetworKit {

using IsomorphismTest::graphOf;
using IsomorphismTest::Match;
using IsomorphismTest::referenceMatches;
using IsomorphismTest::Semantics;

class SubgraphIsomorphismGTest : public testing::Test {
protected:
    /// The interrupt tests set a *process-global* flag. Every one of them clears it again, but an
    /// ASSERT_* returns from the test early, so a future edit could leave it set and silently
    /// interrupt every later test in the binary. Clearing it here makes that impossible.
    void TearDown() override { GlobalState::setReceivedSIGINT(false); }
};

namespace {

/// K_n, which is the cheapest way to get a lot of matches out of a small graph - and a lot of
/// matches is what gives the reporting path a chance to actually collide.
Graph completeGraph(count n) {
    Graph G(n);
    for (node u = 0; u < n; ++u)
        for (node v = u + 1; v < n; ++v)
            G.addEdge(u, v);
    return G;
}

/// Reports a precomputed match set from several threads at once, accumulating per worker exactly
/// as ParallelRI does. Exists so that the base class's parallel reporting path is exercised by
/// something that reports a *known* match set, which keeps a failure here about the base class
/// rather than about whichever search produced the matches.
///
/// This is also the shortest worked example of the parallel half of the protocol: the base class
/// does not collect for you, so a worker buffers its own matches, counts its own, and the merged
/// pair goes to finishRun() once the region has joined.
class MultiWorkerReporter final : public SubgraphIsomorphism {

public:
    MultiWorkerReporter(const Graph &pattern, const Graph &target, Semantics semantics,
                        count numWorkers)
        : SubgraphIsomorphism(pattern, target, semantics, 0), numWorkers(numWorkers) {}

    /// Asked once by run(), and by anybody sizing a per-worker accumulator, so the two cannot
    /// disagree about how many workers there are.
    count numberOfWorkers() const override { return numWorkers; }

    void run() override {
        const std::vector<Match> all = IsomorphismTest::referenceMatches(
            *pattern, *target, semantics, patternNodeLabels, targetNodeLabels);

        Aux::SignalHandler handler;
        prepareRun();

        /// One per worker. Padded so two workers never write to the same cache line, which is
        /// what makes the unsynchronized accumulation below legitimate.
        struct alignas(64) Slot {
            std::vector<Match> buffer;
            count found = 0;
        };
        std::vector<Slot> slots(numWorkers);
        const bool store = storesMatches();

#pragma omp parallel num_threads(static_cast<int>(numWorkers))
        {
            const index tid = static_cast<index>(omp_get_thread_num());

            // The stride is the team OpenMP actually gave, not the one that was asked for.
            // num_threads is an upper bound, so a runtime that hands back fewer threads than
            // requested is legal - and striding by the request would then step over the matches
            // belonging to the threads that never started, losing them silently. `slots` is sized
            // by the request, which is at least the team, so indexing it by tid stays safe.
            // ParallelRIImpl makes the same distinction, under the name activeWorkers.
            const count team = static_cast<count>(omp_get_num_threads());

            Slot &slot = slots[tid];
            // Deal the matches round-robin so every worker reports, and reports interleaved.
            for (index i = tid; i < all.size(); i += team) {
                // Non-throwing inside the region; the throw happens once after the join.
                if (!handler.isRunning())
                    break;
                ++slot.found;
                if (!invokeCallback(tid, all[i]) && store)
                    slot.buffer.push_back(all[i]);
            }
        }

        handler.assureRunning();

        std::vector<Match> merged;
        count found = 0;
        for (Slot &slot : slots) {
            found += slot.found;
            for (Match &match : slot.buffer)
                merged.push_back(std::move(match));
        }

        finishRun(std::move(merged), found);
    }

private:
    count numWorkers;
};

} // namespace

TEST_F(SubgraphIsomorphismGTest, testHarnessDrivesAnAlgorithmThroughTheRunProtocol) {

    // The three helpers in the utils header are templates, so nothing compiles them until somebody
    // calls them. Instantiating them here means a mistake in the harness surfaces now, rather than
    // on the day somebody first tries to use it against a real algorithm.
    //
    // Agreement with the reference is of course trivial for an adapter that *is* the reference.
    // What is not trivial, and what this actually tests, is the protocol around it: prepareRun,
    // reportMatch, finishRun, the three callback forms, setStoreMatches, and the match cap.
    const auto make = [](const Graph &pattern, const Graph &target, Semantics semantics,
                         count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new IsomorphismTest::ReferenceSubgraphIsomorphism(pattern, target, semantics,
                                                              maxMatches));
    };

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
    IsomorphismTest::expectCallbackFormsAgree(make);
}

TEST_F(SubgraphIsomorphismGTest, testSerialCallbackIsNeverEnteredConcurrently) {

    // The guarantee MatchCallback makes, and the one the module got wrong before invokeCallback()
    // took the lock itself: ParallelRIImpl was handed a "serialize this" flag, but the reporting
    // call lives inside RIImpl several levels down, which had no way to act on it. ParallelRIGTest
    // asserts the same thing through the real search; this one pins it to the base class.
    // A 4-path in K8: 8 * 7 * 6 * 5 = 1680 matches, so the four workers really do pile into the
    // callback rather than taking turns by accident.
    const Graph pattern = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});
    const Graph target = completeGraph(8);

    std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);
    ASSERT_EQ(expected.size(), 1680u);

    // The detector: count how many callback invocations are in flight at once. A real user's
    // callback needs none of this - that is precisely the guarantee. The mutex below belongs to
    // the test, not to the guarantee: without it, a regression corrupts the heap and aborts
    // instead of failing an assertion, which is a much worse thing to hand somebody.
    std::atomic<int> inside{0};
    std::atomic<int> maxObserved{0};
    std::mutex collectedMutex;
    std::vector<Match> collected;

    MultiWorkerReporter algo(pattern, target, Semantics::MONOMORPHISM, 4);
    algo.setCallback([&](const Match &match) {
        const int now = inside.fetch_add(1, std::memory_order_acq_rel) + 1;

        int previousMax = maxObserved.load(std::memory_order_relaxed);
        while (now > previousMax
               && !maxObserved.compare_exchange_weak(previousMax, now, std::memory_order_relaxed))
            ;

        {
            const std::lock_guard<std::mutex> guard(collectedMutex);
            collected.push_back(match);
        }

        inside.fetch_sub(1, std::memory_order_acq_rel);
    });
    algo.run();

    // Verified to catch the real thing: with the lock removed from invokeCallback(), this reports
    // 3 or 4 rather than 1.
    EXPECT_EQ(maxObserved.load(), 1) << "the serial callback form was entered concurrently";

    IsomorphismTest::sortMatches(collected);
    EXPECT_EQ(collected, expected);
    EXPECT_EQ(algo.numberOfMatches(), expected.size());
}

TEST_F(SubgraphIsomorphismGTest, testParallelCallbackReceivesEveryMatchOnce) {

    // The other form: invoked directly by each worker with that worker's id, no queueing. Every
    // match must still arrive exactly once, and every tid must be in range.
    const Graph pattern = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});
    const Graph target = completeGraph(8);
    const count numWorkers = 4;

    std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);

    std::vector<std::vector<Match>> perWorker(numWorkers);

    MultiWorkerReporter algo(pattern, target, Semantics::MONOMORPHISM, numWorkers);
    algo.setCallback([&](index tid, const Match &match) {
        ASSERT_LT(tid, numWorkers);
        perWorker[tid].push_back(match); // no lock needed: each tid owns its slot
    });
    algo.run();

    std::vector<Match> collected;
    for (const std::vector<Match> &slot : perWorker)
        collected.insert(collected.end(), slot.begin(), slot.end());

    IsomorphismTest::sortMatches(collected);
    EXPECT_EQ(collected, expected);
    EXPECT_EQ(algo.numberOfMatches(), expected.size());
}

TEST_F(SubgraphIsomorphismGTest, testWorkerBuffersMergeIntoOneResult) {

    // With no callback at all, each worker fills its own padded buffer and the merged pair goes to
    // finishRun(). The result must not depend on how many workers produced it.
    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = graphOf(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});

    std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);

    for (count workers : {count{1}, count{2}, count{4}, count{8}}) {
        MultiWorkerReporter algo(pattern, target, Semantics::MONOMORPHISM, workers);
        algo.run();

        std::vector<Match> actual = algo.getMatches();
        IsomorphismTest::sortMatches(actual);

        EXPECT_EQ(actual, expected) << "workers: " << workers;
        EXPECT_EQ(algo.numberOfMatches(), expected.size()) << "workers: " << workers;
    }
}

TEST_F(SubgraphIsomorphismGTest, testCountingOnlyStoresNothing) {

    // setStoreMatches(false) is what storesMatches() reports to a parallel search, and the whole
    // point is that it then allocates nothing per match. The count must survive regardless.
    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = completeGraph(6);

    const count expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM).size();

    for (count workers : {count{1}, count{4}}) {
        MultiWorkerReporter algo(pattern, target, Semantics::MONOMORPHISM, workers);
        algo.setStoreMatches(false);
        algo.run();

        EXPECT_EQ(algo.numberOfMatches(), expected) << "workers: " << workers;
        EXPECT_TRUE(algo.hasMatch()) << "workers: " << workers;
        EXPECT_THROW(algo.getMatches(), std::runtime_error) << "workers: " << workers;
    }
}

TEST_F(SubgraphIsomorphismGTest, testNumberOfWorkers) {

    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = completeGraph(5);

    // A sequential algorithm reports one worker, before and after the run alike.
    IsomorphismTest::ReferenceSubgraphIsomorphism sequential(pattern, target,
                                                             Semantics::MONOMORPHISM);
    EXPECT_EQ(sequential.numberOfWorkers(), 1u);
    sequential.run();
    EXPECT_EQ(sequential.numberOfWorkers(), 1u);

    // A parallel one reports what it declared - and the answer is available *before* run(), which
    // is the point: it is what a caller sizes a per-worker accumulator by.
    for (count workers : {count{1}, count{2}, count{4}, count{7}}) {
        MultiWorkerReporter algo(pattern, target, Semantics::MONOMORPHISM, workers);
        EXPECT_EQ(algo.numberOfWorkers(), workers);
        algo.run();
        EXPECT_EQ(algo.numberOfWorkers(), workers);
    }
}

TEST_F(SubgraphIsomorphismGTest, testInterruptLeavesTheAlgorithmUnfinished) {

    // The documented policy: run() throws, the algorithm is left not-finished, and every result
    // accessor throws. Driven by setting the global flag rather than by raising a real SIGINT,
    // which would be untestable.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = completeGraph(6);

    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    count delivered = 0;
    algo.setCallback([&](const Match &) {
        // Interrupt partway through, so the search is genuinely mid-flight.
        if (++delivered == 5)
            GlobalState::setReceivedSIGINT(true);
    });

    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);

    EXPECT_FALSE(algo.hasFinished()) << "an interrupted run must not count as finished";
    EXPECT_THROW(algo.numberOfMatches(), std::runtime_error);
    EXPECT_THROW(algo.hasMatch(), std::runtime_error);

    // Matches already handed to the callback stay handed over - a search cannot take them back.
    EXPECT_GE(delivered, 5u);
}

TEST_F(SubgraphIsomorphismGTest, testAlgorithmRecoversAfterAnInterrupt) {

    // An interrupted object must not be poisoned: a second, clean run has to give the full answer.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = completeGraph(6);

    std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);

    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    GlobalState::setReceivedSIGINT(true);
    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);
    ASSERT_FALSE(algo.hasFinished());

    algo.run();

    ASSERT_TRUE(algo.hasFinished());
    std::vector<Match> actual = algo.getMatches();
    IsomorphismTest::sortMatches(actual);
    EXPECT_EQ(actual, expected);
}

TEST_F(SubgraphIsomorphismGTest, testSetEdgeLabelsValidatesItsInput) {

    // The setter is the whole check: nothing consumes the vectors until run(), so a vector that is
    // too short is a silent out-of-bounds read in the innermost loop of whichever search gets it.
    Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});

    IsomorphismTest::ReferenceSubgraphIsomorphism unindexed(pattern, target,
                                                            Semantics::MONOMORPHISM);
    ASSERT_FALSE(pattern.hasEdgeIds());
    EXPECT_THROW(unindexed.setEdgeLabels({1, 2}, {1, 2, 3}), std::runtime_error)
        << "edge labels are indexed by edge id, so a graph without ids has no index space";

    // Clearing must keep working even then: it is tested before the size checks, which is the only
    // reason the documented "pass two empty vectors" idiom survives them.
    EXPECT_NO_THROW(unindexed.setEdgeLabels({}, {}));

    pattern.indexEdges();
    target.indexEdges();
    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    EXPECT_THROW(algo.setEdgeLabels({1}, {1, 2, 3}), std::runtime_error) << "pattern too short";
    EXPECT_THROW(algo.setEdgeLabels({1, 2}, {1, 2}), std::runtime_error) << "target too short";
    EXPECT_NO_THROW(algo.setEdgeLabels({1, 2}, {1, 2, 3}));
    EXPECT_NO_THROW(algo.setEdgeLabels({}, {}));
}

TEST_F(SubgraphIsomorphismGTest, testEdgeLabelsAreEitherHonouredOrRefused) {

    // The module's rule, judged here rather than in any one algorithm's own test file because it
    // applies to every driver: an algorithm either honours the edge labels it was given or refuses
    // them outright. Silently returning matches that violate an edge label the caller asked for is
    // the one failure worse than refusing.
    const IsomorphismTest::LabelledGraph pattern =
        IsomorphismTest::labelledGraphOf(3, {{0, 1, 1}, {1, 2, 2}});
    const IsomorphismTest::LabelledGraph target =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}});

    // The refusing side is VF3 alone. Its search is still unwritten, so honouring a label is not
    // something it could do even if it tried.
    VF3 vf3(pattern.G, target.G, Semantics::MONOMORPHISM);
    vf3.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    EXPECT_THROW(vf3.run(), std::runtime_error);
    EXPECT_FALSE(vf3.hasFinished()) << "a refused run must not count as finished";

    // The honouring side is VF2, RI and ParallelRI. Each must answer the labelled question, not
    // the unlabelled one, which is why every count is compared against the reference rather than
    // merely being nonzero.
    const count labelled = referenceMatches(pattern.G, target.G, Semantics::MONOMORPHISM, {}, {},
                                            pattern.edgeLabels, target.edgeLabels)
                               .size();
    const count unlabelled = referenceMatches(pattern.G, target.G, Semantics::MONOMORPHISM).size();
    ASSERT_LT(labelled, unlabelled) << "the labels must rule some match out, or the comparison "
                                       "below cannot tell an honoured label from an ignored one";

    VF2 vf2(pattern.G, target.G, Semantics::MONOMORPHISM);
    vf2.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    EXPECT_NO_THROW(vf2.run());
    EXPECT_TRUE(vf2.hasFinished());
    EXPECT_EQ(vf2.numberOfMatches(), labelled);

    RI ri(pattern.G, target.G, RI::Variant::RI, Semantics::MONOMORPHISM);
    ri.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    EXPECT_NO_THROW(ri.run());
    EXPECT_TRUE(ri.hasFinished());
    EXPECT_EQ(ri.numberOfMatches(), labelled);

    ParallelRI parallelRi(pattern.G, target.G, RI::Variant::RI, Semantics::MONOMORPHISM);
    parallelRi.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    EXPECT_NO_THROW(parallelRi.run());
    EXPECT_TRUE(parallelRi.hasFinished());
    EXPECT_EQ(parallelRi.numberOfMatches(), labelled);

    // Passing no labels at all must still reach the larger answer, so the labelled count above is
    // a restriction the search applied and not a search that lost matches for another reason.
    VF2 vf2Unlabelled(pattern.G, target.G, Semantics::MONOMORPHISM);
    EXPECT_NO_THROW(vf2Unlabelled.run());
    EXPECT_TRUE(vf2Unlabelled.hasFinished());
    EXPECT_EQ(vf2Unlabelled.numberOfMatches(), unlabelled);
}

TEST_F(SubgraphIsomorphismGTest, testParallelEdgesWithDisagreeingLabelsAreRefused) {

    // The other half of the module's rule: the algorithms that *will* understand edge labels refuse
    // only what no snapshot can represent. The two parallel 0-1 edges disagree, so collapsing them
    // leaves one arc that cannot carry both labels.
    //
    // Both algorithms refuse before any search begins, and ParallelRI refuses before any thread
    // is started, so there is never a worker left to unwind. The equally-labelled case below is
    // what proves the refusal is this narrow one and not a blanket rejection of edge labels.
    const IsomorphismTest::LabelledGraph pattern =
        IsomorphismTest::labelledGraphOf(3, {{0, 1, 1}, {1, 2, 2}});
    const IsomorphismTest::LabelledGraph target =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {0, 1, 4}, {1, 2, 2}, {2, 3, 1}});

    RI ri(pattern.G, target.G, RI::Variant::RI, Semantics::MONOMORPHISM);
    ri.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    EXPECT_THROW(ri.run(), std::runtime_error);
    EXPECT_FALSE(ri.hasFinished());

    ParallelRI parallelRi(pattern.G, target.G, RI::Variant::RI, Semantics::MONOMORPHISM);
    parallelRi.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    EXPECT_THROW(parallelRi.run(), std::runtime_error);
    EXPECT_FALSE(parallelRi.hasFinished());

    // Equally-labelled parallel edges collapse losslessly, so they must get past the refusal and
    // be searched normally.
    const IsomorphismTest::LabelledGraph agreeing =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {0, 1, 1}, {1, 2, 2}, {2, 3, 1}});

    RI lossless(pattern.G, agreeing.G, RI::Variant::RI, Semantics::MONOMORPHISM);
    lossless.setEdgeLabels(pattern.edgeLabels, agreeing.edgeLabels);
    EXPECT_NO_THROW(lossless.run())
        << "collapsing equally-labelled parallel edges is lossless and must not be refused";
    EXPECT_TRUE(lossless.hasFinished());
    EXPECT_EQ(lossless.numberOfMatches(),
              referenceMatches(pattern.G, agreeing.G, Semantics::MONOMORPHISM, {}, {},
                               pattern.edgeLabels, agreeing.edgeLabels)
                  .size());
}

// ------------------------------------------------------------------------------------------
// Revalidation at run()
//
// The two graphs are held by pointer, so everything checked at construction or configuration
// time can go stale before run() is reached. These pin down that run() rechecks rather than
// trusting what was true earlier - the label cases especially, where a short vector is an
// out-of-bounds read rather than merely a wrong answer.
// ------------------------------------------------------------------------------------------

TEST_F(SubgraphIsomorphismGTest, testRunRechecksGraphInvariants) {

    Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    ASSERT_NO_THROW(algo.run());

    // The constructor rejected self-loops in the pattern. Adding one afterwards must not sneak
    // an undefined query past that check.
    pattern.addEdge(1, 1);
    ASSERT_EQ(pattern.numberOfSelfLoops(), 1u);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "a self-loop added to the pattern after construction must still be rejected";
    EXPECT_FALSE(algo.hasFinished());
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksNodeLabelSizes) {

    const Graph pattern = graphOf(2, {{0, 1}});
    Graph target = graphOf(3, {{0, 1}, {1, 2}});

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    algo.setNodeLabels({1, 1}, {1, 1, 1});
    ASSERT_NO_THROW(algo.run());

    // The label vector was sized against the target as it was. One more node and it is short,
    // which without the recheck is a read past its end rather than an error.
    target.addNode();
    ASSERT_GT(target.upperNodeIdBound(), 3u);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "a node added after setNodeLabels() leaves the label vector short";
    EXPECT_FALSE(algo.hasFinished());

    // Re-supplying the labels at the new size is all it takes to make the run legal again.
    algo.setNodeLabels({1, 1}, {1, 1, 1, 1});
    EXPECT_NO_THROW(algo.run());
    EXPECT_TRUE(algo.hasFinished());
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksEdgeLabelSizes) {

    const IsomorphismTest::LabelledGraph pattern = IsomorphismTest::labelledGraphOf(2, {{0, 1, 1}});
    IsomorphismTest::LabelledGraph target =
        IsomorphismTest::labelledGraphOf(3, {{0, 1, 1}, {1, 2, 1}});

    // RI honours edge labels, so it reaches the size check rather than refusing them outright
    // the way VF2 does. Unlike the node label case above, this one has a second line of defence:
    // SearchGraph's constructor validates the edge label vector while building the snapshot, so
    // it would throw here even without prepareRun(). The test pins the guarantee at the run()
    // boundary rather than to whichever layer happens to enforce it.
    RI algo(pattern.G, target.G, RI::Variant::RI, Semantics::MONOMORPHISM);
    algo.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    ASSERT_NO_THROW(algo.run());

    const index boundBefore = target.G.upperEdgeIdBound();
    target.G.addEdge(0, 2);
    target.G.indexEdges(true);
    ASSERT_GT(target.G.upperEdgeIdBound(), boundBefore);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "an edge added after setEdgeLabels() leaves the edge label vector short";
    EXPECT_FALSE(algo.hasFinished());

    std::vector<index> grown = target.edgeLabels;
    grown.resize(target.G.upperEdgeIdBound(), 1);
    algo.setEdgeLabels(pattern.edgeLabels, grown);
    EXPECT_NO_THROW(algo.run());
    EXPECT_TRUE(algo.hasFinished());
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksLeaveNothingQueryableWhenTheyFire) {

    const Graph pattern = graphOf(2, {{0, 1}});
    Graph target = graphOf(3, {{0, 1}, {1, 2}});

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    algo.setNodeLabels({1, 1}, {1, 1, 1});
    ASSERT_NO_THROW(algo.run());
    ASSERT_GT(algo.numberOfMatches(), 0u);

    target.addNode();
    ASSERT_THROW(algo.run(), std::runtime_error);

    // The recheck runs after the reset on purpose: a rejected run must not hand back the matches
    // of the run before it.
    EXPECT_FALSE(algo.hasFinished());
    EXPECT_THROW(algo.getMatches(), std::runtime_error);
}

} // namespace NetworKit
