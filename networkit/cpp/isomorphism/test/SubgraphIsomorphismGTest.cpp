/*
 * SubgraphIsomorphismGTest.cpp
 *
 * Tests the base class SubgraphIsomorphism: the run() protocol, the two callback forms, the
 * per-worker accumulation of a parallel search, input validation and interruption. Most tests use
 * two stand-ins that report a precomputed match set: ReferenceSubgraphIsomorphism for the
 * sequential and MultiWorkerReporter for the parallel path.
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

#include "SubgraphIsomorphismTestUtils.hpp"

namespace NetworKit {

using IsomorphismTest::graphOf;
using IsomorphismTest::Match;
using IsomorphismTest::referenceMatches;
using IsomorphismTest::Semantics;

class SubgraphIsomorphismGTest : public testing::Test {
protected:
    /// The interrupt tests set a process-global flag, which a failing ASSERT_* would leave set.
    void TearDown() override { GlobalState::setReceivedSIGINT(false); }
};

namespace {

Graph completeGraph(count n) {
    Graph G(n);
    for (node u = 0; u < n; ++u)
        for (node v = u + 1; v < n; ++v)
            G.addEdge(u, v);
    return G;
}

/// Reports the reference matches from several threads, as ParallelRI does. It ignores
/// @a maxMatches, like workers that overshoot the cap.
class MultiWorkerReporter final : public SubgraphIsomorphism {

public:
    MultiWorkerReporter(const Graph &pattern, const Graph &target, Semantics semantics,
                        count numWorkers, count maxMatches = 0)
        : SubgraphIsomorphism(pattern, target, semantics, maxMatches), numWorkers(numWorkers) {}

    count numberOfWorkers() const override { return numWorkers; }

    void run() override {
        const std::vector<Match> all = IsomorphismTest::referenceMatches(
            *pattern, *target, semantics, patternNodeLabels, targetNodeLabels);

        Aux::SignalHandler handler;
        prepareRun();

        struct alignas(64) Slot {
            std::vector<Match> buffer;
            count found = 0;
        };
        std::vector<Slot> slots(numWorkers);
        const bool store = storesMatches();

#pragma omp parallel num_threads(static_cast<int>(numWorkers))
        {
            const index tid = static_cast<index>(omp_get_thread_num());

            // Stride by the team size, since OpenMP may start fewer threads than requested.
            const count team = static_cast<count>(omp_get_num_threads());

            Slot &slot = slots[tid];
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

TEST_F(SubgraphIsomorphismGTest, testBothCallbackFormsReceiveEveryMatchOnce) {

    // A 4-path in K8 has 1680 matches, so the four workers compete for the callback.
    const Graph pattern = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});
    const Graph target = completeGraph(8);
    const count numWorkers = 4;

    std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);
    ASSERT_EQ(expected.size(), 1680u);

    // The mutex keeps a regression from corrupting `collected` instead of failing the assertion.
    std::atomic<int> inside{0};
    std::atomic<int> maxObserved{0};
    std::mutex collectedMutex;
    std::vector<Match> collected;

    MultiWorkerReporter serial(pattern, target, Semantics::MONOMORPHISM, numWorkers);
    serial.setCallback([&](const Match &match) {
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
    serial.run();

    EXPECT_EQ(maxObserved.load(), 1) << "the serial callback form was entered concurrently";

    IsomorphismTest::sortMatches(collected);
    EXPECT_EQ(collected, expected);
    EXPECT_EQ(serial.numberOfMatches(), expected.size());

    std::vector<std::vector<Match>> perWorker(numWorkers);

    MultiWorkerReporter parallel(pattern, target, Semantics::MONOMORPHISM, numWorkers);
    parallel.setCallback([&](index tid, const Match &match) {
        ASSERT_LT(tid, numWorkers);
        perWorker[tid].push_back(match); // no lock needed: each tid owns its slot
    });
    parallel.run();

    std::vector<Match> fromWorkers;
    for (const std::vector<Match> &slot : perWorker)
        fromWorkers.insert(fromWorkers.end(), slot.begin(), slot.end());

    IsomorphismTest::sortMatches(fromWorkers);
    EXPECT_EQ(fromWorkers, expected);
    EXPECT_EQ(parallel.numberOfMatches(), expected.size());
}

TEST_F(SubgraphIsomorphismGTest, testTheAnswerDoesNotDependOnTheWorkerCount) {

    const Graph pattern = graphOf(2, {{0, 1}});

    const Graph cycle = graphOf(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});

    std::vector<Match> expected = referenceMatches(pattern, cycle, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);

    for (count workers : {count{1}, count{2}, count{4}, count{8}}) {
        MultiWorkerReporter algo(pattern, cycle, Semantics::MONOMORPHISM, workers);
        algo.run();

        std::vector<Match> actual = algo.getMatches();
        IsomorphismTest::sortMatches(actual);

        EXPECT_EQ(actual, expected) << "workers: " << workers;
        EXPECT_EQ(algo.numberOfMatches(), expected.size()) << "workers: " << workers;
    }

    const Graph k6 = completeGraph(6);
    const count expectedCount = referenceMatches(pattern, k6, Semantics::MONOMORPHISM).size();

    for (count workers : {count{1}, count{4}}) {
        MultiWorkerReporter algo(pattern, k6, Semantics::MONOMORPHISM, workers);
        algo.setStoreMatches(false);
        algo.run();

        EXPECT_EQ(algo.numberOfMatches(), expectedCount) << "workers: " << workers;
        EXPECT_TRUE(algo.hasMatch()) << "workers: " << workers;
        EXPECT_THROW(algo.getMatches(), std::runtime_error) << "workers: " << workers;
    }

    const Graph k5 = completeGraph(5);

    IsomorphismTest::ReferenceSubgraphIsomorphism sequential(pattern, k5, Semantics::MONOMORPHISM);
    EXPECT_EQ(sequential.numberOfWorkers(), 1u);
    sequential.run();
    EXPECT_EQ(sequential.numberOfWorkers(), 1u);

    for (count workers : {count{1}, count{2}, count{4}, count{7}}) {
        MultiWorkerReporter algo(pattern, k5, Semantics::MONOMORPHISM, workers);
        EXPECT_EQ(algo.numberOfWorkers(), workers);
        algo.run();
        EXPECT_EQ(algo.numberOfWorkers(), workers);
    }
}

TEST_F(SubgraphIsomorphismGTest, testOvershootIsTrimmedUnlessACallbackReceivedIt) {

    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = completeGraph(5);
    const count total = referenceMatches(pattern, target, Semantics::MONOMORPHISM).size();
    constexpr count cap = 5;
    constexpr count numWorkers = 4;
    ASSERT_GT(total, cap);

    MultiWorkerReporter stored(pattern, target, Semantics::MONOMORPHISM, numWorkers, cap);
    stored.run();
    EXPECT_EQ(stored.numberOfMatches(), cap);
    EXPECT_EQ(stored.getMatches().size(), cap);

    MultiWorkerReporter counted(pattern, target, Semantics::MONOMORPHISM, numWorkers, cap);
    counted.setStoreMatches(false);
    counted.run();
    EXPECT_EQ(counted.numberOfMatches(), cap);

    count delivered = 0;
    MultiWorkerReporter called(pattern, target, Semantics::MONOMORPHISM, numWorkers, cap);
    called.setCallback([&](const Match &) { ++delivered; });
    called.run();
    EXPECT_EQ(delivered, total);
    EXPECT_EQ(called.numberOfMatches(), total) << "delivered matches cannot be taken back";
    EXPECT_THROW(called.getMatches(), std::runtime_error);
}

TEST_F(SubgraphIsomorphismGTest, testInterruptLeavesTheAlgorithmUnfinishedButUsable) {

    // The test sets the global flag instead of raising a real SIGINT.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = completeGraph(6);

    std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);

    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    count delivered = 0;
    algo.setCallback([&](const Match &) {
        if (++delivered == 5)
            GlobalState::setReceivedSIGINT(true);
    });

    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);

    EXPECT_FALSE(algo.hasFinished()) << "an interrupted run must not count as finished";
    EXPECT_THROW(algo.numberOfMatches(), std::runtime_error);
    EXPECT_THROW(algo.hasMatch(), std::runtime_error);

    EXPECT_GE(delivered, 5u);

    // A second object without a callback, since the recovery check needs the stored matches.
    IsomorphismTest::ReferenceSubgraphIsomorphism recovering(pattern, target,
                                                             Semantics::MONOMORPHISM);

    GlobalState::setReceivedSIGINT(true);
    EXPECT_THROW(recovering.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);
    ASSERT_FALSE(recovering.hasFinished());

    recovering.run();

    ASSERT_TRUE(recovering.hasFinished());
    std::vector<Match> actual = recovering.getMatches();
    IsomorphismTest::sortMatches(actual);
    EXPECT_EQ(actual, expected);
}

TEST_F(SubgraphIsomorphismGTest, testGraphsMustAgreeOnDirectedness) {

    const Graph undirected = graphOf(2, {{0, 1}});
    const Graph directed = graphOf(2, {{0, 1}}, true);

    EXPECT_THROW(VF2 algo(undirected, directed), std::runtime_error);
    EXPECT_THROW(VF2 algo(directed, undirected), std::runtime_error);
    EXPECT_NO_THROW(VF2 algo(directed, directed));
}

TEST_F(SubgraphIsomorphismGTest, testSetNodeLabelsValidatesItsInput) {

    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = graphOf(3, {{0, 1}, {1, 2}});
    VF2 algo(pattern, target, Semantics::MONOMORPHISM);

    EXPECT_THROW(algo.setNodeLabels({1}, {1, 1, 1}), std::runtime_error) << "pattern too short";
    EXPECT_THROW(algo.setNodeLabels({1, 1}, {1, 1}), std::runtime_error) << "target too short";
    EXPECT_NO_THROW(algo.setNodeLabels({1, 1}, {1, 1, 1}));
    EXPECT_NO_THROW(algo.setNodeLabels({}, {}));
}

TEST_F(SubgraphIsomorphismGTest, testSetEdgeLabelsValidatesItsInput) {

    Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});

    IsomorphismTest::ReferenceSubgraphIsomorphism unindexed(pattern, target,
                                                            Semantics::MONOMORPHISM);
    ASSERT_FALSE(pattern.hasEdgeIds());
    EXPECT_THROW(unindexed.setEdgeLabels({1, 2}, {1, 2, 3}), std::runtime_error)
        << "edge labels are indexed by edge id, so a graph without ids has no index space";

    // Two empty vectors clear the labels, even without edge ids.
    EXPECT_NO_THROW(unindexed.setEdgeLabels({}, {}));

    pattern.indexEdges();
    EXPECT_THROW(unindexed.setEdgeLabels({1, 2}, {1, 2, 3}), std::runtime_error)
        << "the target needs edge ids as well";

    target.indexEdges();
    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    EXPECT_THROW(algo.setEdgeLabels({1}, {1, 2, 3}), std::runtime_error) << "pattern too short";
    EXPECT_THROW(algo.setEdgeLabels({1, 2}, {1, 2}), std::runtime_error) << "target too short";
    EXPECT_NO_THROW(algo.setEdgeLabels({1, 2}, {1, 2, 3}));
    EXPECT_NO_THROW(algo.setEdgeLabels({}, {}));
}

TEST_F(SubgraphIsomorphismGTest, testEdgeLabelsAreHonoured) {

    const IsomorphismTest::LabelledGraph pattern =
        IsomorphismTest::labelledGraphOf(3, {{0, 1, 1}, {1, 2, 2}});
    const IsomorphismTest::LabelledGraph target =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}});

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

    VF2 vf2Unlabelled(pattern.G, target.G, Semantics::MONOMORPHISM);
    EXPECT_NO_THROW(vf2Unlabelled.run());
    EXPECT_TRUE(vf2Unlabelled.hasFinished());
    EXPECT_EQ(vf2Unlabelled.numberOfMatches(), unlabelled);
}

TEST_F(SubgraphIsomorphismGTest, testParallelEdgesWithDisagreeingLabelsAreRefused) {

    // One arc of the snapshot cannot represent the two differently labelled 0-1 edges.
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

TEST_F(SubgraphIsomorphismGTest, testRunRechecksGraphInvariants) {

    Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    ASSERT_NO_THROW(algo.run());

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
    ASSERT_GT(algo.numberOfMatches(), 0u);

    target.addNode();
    ASSERT_GT(target.upperNodeIdBound(), 3u);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "a node added after setNodeLabels() leaves the label vector short";

    // A rejected run must not hand back the matches of the previous run.
    EXPECT_FALSE(algo.hasFinished());
    EXPECT_THROW(algo.getMatches(), std::runtime_error);

    algo.setNodeLabels({1, 1}, {1, 1, 1, 1});
    EXPECT_NO_THROW(algo.run());
    EXPECT_TRUE(algo.hasFinished());
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksEdgeLabelSizes) {

    const IsomorphismTest::LabelledGraph pattern = IsomorphismTest::labelledGraphOf(2, {{0, 1, 1}});
    IsomorphismTest::LabelledGraph target =
        IsomorphismTest::labelledGraphOf(3, {{0, 1, 1}, {1, 2, 1}});

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

} // namespace NetworKit
