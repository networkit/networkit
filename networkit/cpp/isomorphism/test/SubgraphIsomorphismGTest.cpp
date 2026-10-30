/*
 * SubgraphIsomorphismGTest.cpp
 *
 * Tests the base class SubgraphIsomorphism: the run() protocol, the two callback forms, the
 * per-worker accumulation of a parallel search, input validation and interruption. Most tests use
 * two stand-ins that report a precomputed match set: ReferenceSubgraphIsomorphism for the
 * sequential and MultiWorkerReporter for the parallel path.
 */

#include <algorithm>
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <omp.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/GlobalState.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/ParallelRI.hpp>
#include <networkit/isomorphism/RI.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>
#include <networkit/isomorphism/VF2.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"

namespace NetworKit {
namespace {

using IsomorphismTest::completeGraph;
using IsomorphismTest::graphOf;
using IsomorphismTest::Match;
using IsomorphismTest::referenceMatches;
using IsomorphismTest::Semantics;
using testing::HasSubstr;
using testing::SizeIs;
using testing::ThrowsMessage;
using testing::UnorderedElementsAreArray;

class SubgraphIsomorphismGTest : public testing::Test {
protected:
    /// The interrupt tests set a process-global flag, which a failing ASSERT_* would leave set.
    void TearDown() override { GlobalState::setReceivedSIGINT(false); }
};

/// Reports the reference matches from several threads, as a parallel algorithm does. It ignores
/// @a maxMatches, like workers that overshoot the cap.
class MultiWorkerReporter final : public SubgraphIsomorphism {

public:
    MultiWorkerReporter(const Graph &pattern, const Graph &target, Semantics semantics,
                        count numWorkers, count maxMatches = 0)
        : SubgraphIsomorphism(pattern, target, semantics, maxMatches), numWorkers(numWorkers) {}

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

        // libomp leaks memory when a team grows past twice the CPU count, which fails the leak
        // checks. We therefore request at most that many threads.
        const count threads = std::min(numWorkers, 2 * static_cast<count>(omp_get_num_procs()));

#pragma omp parallel num_threads(static_cast<int>(threads))
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

using Factory = std::function<std::unique_ptr<SubgraphIsomorphism>(const Graph &, const Graph &)>;

/// The real algorithms, for the base-class behavior that each of them has to keep.
std::vector<std::pair<std::string, Factory>> algorithms() {
    return {{"VF2",
             [](const Graph &pattern, const Graph &target) -> std::unique_ptr<SubgraphIsomorphism> {
                 return std::make_unique<VF2>(pattern, target, Semantics::MONOMORPHISM);
             }},
            {"RI",
             [](const Graph &pattern, const Graph &target) -> std::unique_ptr<SubgraphIsomorphism> {
                 return std::make_unique<RI>(pattern, target, RI::Variant::RI,
                                             Semantics::MONOMORPHISM);
             }},
            {"ParallelRI",
             [](const Graph &pattern, const Graph &target) -> std::unique_ptr<SubgraphIsomorphism> {
                 return std::make_unique<ParallelRI>(pattern, target, RI::Variant::RI,
                                                     Semantics::MONOMORPHISM);
             }}};
}

TEST_F(SubgraphIsomorphismGTest, testHarnessDrivesAnAlgorithmThroughTheRunProtocol) {
    const auto make = [](const Graph &pattern, const Graph &target, Semantics semantics,
                         count maxMatches) -> std::unique_ptr<SubgraphIsomorphism> {
        return std::make_unique<IsomorphismTest::ReferenceSubgraphIsomorphism>(
            pattern, target, semantics, maxMatches);
    };

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
}

TEST_F(SubgraphIsomorphismGTest, testBothCallbackFormsReceiveEveryMatchOnce) {
    // A 4-path in K8 has 1680 matches, so the four workers compete for the callback.
    const Graph pattern = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});
    const Graph target = completeGraph(8);
    const count numWorkers = 4;

    const std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    ASSERT_THAT(expected, SizeIs(1680));

    // The mutex keeps a regression from corrupting `collected` instead of failing the assertion.
    std::atomic<int> inside{0};
    std::atomic<int> maxObserved{0};
    std::mutex collectedMutex;
    std::vector<Match> collected;

    MultiWorkerReporter serial(pattern, target, Semantics::MONOMORPHISM, numWorkers);
    serial.setCallback([&](std::span<const node> match) {
        const int now = inside.fetch_add(1, std::memory_order_acq_rel) + 1;
        Aux::Parallel::atomic_max(maxObserved, now);

        {
            const std::lock_guard<std::mutex> guard(collectedMutex);
            collected.emplace_back(match.begin(), match.end());
        }

        inside.fetch_sub(1, std::memory_order_acq_rel);
    });
    serial.run();

    EXPECT_EQ(maxObserved.load(), 1) << "the serial callback form was entered concurrently";

    EXPECT_THAT(collected, UnorderedElementsAreArray(expected));
    EXPECT_EQ(serial.numberOfMatches(), expected.size());

    std::vector<std::vector<Match>> perWorker(numWorkers);

    MultiWorkerReporter parallel(pattern, target, Semantics::MONOMORPHISM, numWorkers);
    parallel.setParallelCallback([&](index tid, std::span<const node> match) {
        ASSERT_LT(tid, numWorkers);
        perWorker[tid].emplace_back(match.begin(), match.end()); // each tid owns its slot
    });
    parallel.run();

    std::vector<Match> fromWorkers;
    for (const std::vector<Match> &slot : perWorker)
        fromWorkers.insert(fromWorkers.end(), slot.begin(), slot.end());

    EXPECT_THAT(fromWorkers, UnorderedElementsAreArray(expected));
    EXPECT_EQ(parallel.numberOfMatches(), expected.size());
}

TEST_F(SubgraphIsomorphismGTest, testTheAnswerDoesNotDependOnTheWorkerCount) {
    const Graph pattern = graphOf(2, {{0, 1}});

    const Graph cycle = graphOf(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});

    const std::vector<Match> expected = referenceMatches(pattern, cycle, Semantics::MONOMORPHISM);

    for (count workers : {count{1}, count{2}, count{4}, count{8}}) {
        MultiWorkerReporter algo(pattern, cycle, Semantics::MONOMORPHISM, workers);
        algo.run();

        EXPECT_THAT(algo.getMatches(), UnorderedElementsAreArray(expected))
            << "workers: " << workers;
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
    EXPECT_THAT(stored.getMatches(), SizeIs(cap));

    MultiWorkerReporter counted(pattern, target, Semantics::MONOMORPHISM, numWorkers, cap);
    counted.setStoreMatches(false);
    counted.run();
    EXPECT_EQ(counted.numberOfMatches(), cap);

    count delivered = 0;
    MultiWorkerReporter called(pattern, target, Semantics::MONOMORPHISM, numWorkers, cap);
    called.setCallback([&](std::span<const node>) { ++delivered; });
    called.run();
    EXPECT_EQ(delivered, total);
    EXPECT_EQ(called.numberOfMatches(), total) << "delivered matches cannot be taken back";
    EXPECT_THROW(called.getMatches(), std::runtime_error);
}

TEST_F(SubgraphIsomorphismGTest, testInterruptLeavesTheAlgorithmUnfinishedButUsable) {
    // The test sets the global flag instead of raising a real SIGINT.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = completeGraph(6);

    const std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);

    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    count delivered = 0;
    algo.setCallback([&](std::span<const node>) {
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
    EXPECT_THAT(recovering.getMatches(), UnorderedElementsAreArray(expected));
}

TEST_F(SubgraphIsomorphismGTest, testGraphsMustAgreeOnDirectedness) {
    const Graph undirected = graphOf(2, {{0, 1}});
    const Graph directed = graphOf(2, {{0, 1}}, true);

    // The graphs may change after construction, so only run() checks them.
    VF2 undirectedPattern(undirected, directed);
    EXPECT_THROW(undirectedPattern.run(), std::runtime_error);

    VF2 directedPattern(directed, undirected);
    EXPECT_THROW(directedPattern.run(), std::runtime_error);

    VF2 agreeing(directed, directed);
    EXPECT_NO_THROW(agreeing.run());
}

TEST_F(SubgraphIsomorphismGTest, testUnknownSemanticsIsRejected) {
    // Python passes the semantics as a plain integer, so any value can arrive here.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    EXPECT_THROW(VF2 algo(pattern, target, static_cast<Semantics>(2)), std::invalid_argument);
    EXPECT_THROW(VF2 algo(pattern, target, static_cast<Semantics>(256)), std::invalid_argument);
    EXPECT_NO_THROW(VF2 algo(pattern, target, Semantics::MONOMORPHISM));
}

/// ParallelRI hands its matches over only after the join, so whether a run stores them has to be
/// fixed by prepareRun() rather than read from the settings afterwards.
TEST_F(SubgraphIsomorphismGTest, testGetMatchesFollowsTheLastRun) {
    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});

    for (const auto &[name, make] : algorithms()) {
        SCOPED_TRACE(name);

        const std::unique_ptr<SubgraphIsomorphism> algo = make(pattern, target);
        algo->setStoreMatches(false);
        algo->run();
        ASSERT_EQ(algo->numberOfMatches(), 8u);

        // The run only counted, so turning storage on afterwards must not expose an empty result.
        algo->setStoreMatches(true);
        EXPECT_THROW(algo->getMatches(), std::runtime_error);

        algo->run();
        EXPECT_THAT(algo->getMatches(), SizeIs(8));

        // A callback set after the run does not take away the matches that run stored.
        algo->setCallback([](std::span<const node>) {});
        EXPECT_THAT(algo->getMatches(), SizeIs(8));
    }
}

TEST_F(SubgraphIsomorphismGTest, testRunValidatesNodeLabels) {
    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = graphOf(3, {{0, 1}, {1, 2}});
    const std::vector<index> one{1}, two{1, 1}, three{1, 1, 1};

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);

    algo.setNodeLabels(one, three);
    EXPECT_THROW(algo.run(), std::runtime_error) << "pattern too short";

    algo.setNodeLabels(two, two);
    EXPECT_THROW(algo.run(), std::runtime_error) << "target too short";

    algo.setNodeLabels(two, three);
    EXPECT_NO_THROW(algo.run());

    algo.setNodeLabels({}, {});
    EXPECT_NO_THROW(algo.run());
}

TEST_F(SubgraphIsomorphismGTest, testRunValidatesEdgeLabels) {
    Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});
    const std::vector<index> one{1}, two{1, 2}, three{1, 2, 3};

    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    // Two empty vectors clear the labels, even without edge ids.
    algo.setEdgeLabels({}, {});
    EXPECT_NO_THROW(algo.run());

    algo.setEdgeLabels(two, three);
    ASSERT_FALSE(pattern.hasEdgeIds());
    EXPECT_THROW(algo.run(), std::runtime_error)
        << "edge labels are indexed by edge id, so a graph without ids has no index space";

    pattern.indexEdges();
    EXPECT_THROW(algo.run(), std::runtime_error) << "the target needs edge ids as well";

    target.indexEdges();
    EXPECT_NO_THROW(algo.run());

    algo.setEdgeLabels(one, three);
    EXPECT_THROW(algo.run(), std::runtime_error) << "pattern too short";

    algo.setEdgeLabels(two, two);
    EXPECT_THROW(algo.run(), std::runtime_error) << "target too short";
}

TEST_F(SubgraphIsomorphismGTest, testParallelEdgesWithDisagreeingLabelsAreRefused) {
    // One arc of the snapshot cannot represent the two differently labelled 0-1 edges.
    const IsomorphismTest::LabelledGraph plain =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}});
    const IsomorphismTest::LabelledGraph disagreeing =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {0, 1, 4}, {1, 2, 2}, {2, 3, 1}});

    for (const auto &[name, make] : algorithms()) {
        SCOPED_TRACE(name);
        for (const std::string_view side : {"pattern", "target"}) {
            const IsomorphismTest::LabelledGraph &pattern = side == "pattern" ? disagreeing : plain;
            const IsomorphismTest::LabelledGraph &target = side == "pattern" ? plain : disagreeing;

            std::unique_ptr<SubgraphIsomorphism> algo = make(pattern.G, target.G);
            algo->setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
            EXPECT_THAT([&] { algo->run(); }, ThrowsMessage<std::runtime_error>(HasSubstr(side)));
            EXPECT_FALSE(algo->hasFinished());
        }
    }
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

    const std::vector<index> patternLabels{1, 1};

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    algo.setNodeLabels(patternLabels, std::vector<index>{1, 1, 1});
    ASSERT_NO_THROW(algo.run());
    ASSERT_GT(algo.numberOfMatches(), 0u);

    target.addNode();
    ASSERT_GT(target.upperNodeIdBound(), 3u);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "a node added after setNodeLabels() leaves the label vector short";

    // A rejected run must not hand back the matches of the previous run.
    EXPECT_FALSE(algo.hasFinished());
    EXPECT_THROW(algo.getMatches(), std::runtime_error);

    algo.setNodeLabels(patternLabels, std::vector<index>{1, 1, 1, 1});
    EXPECT_NO_THROW(algo.run());
    EXPECT_TRUE(algo.hasFinished());
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksEdgeLabelSizes) {
    const IsomorphismTest::LabelledGraph pattern = IsomorphismTest::labelledGraphOf(2, {{0, 1, 1}});
    IsomorphismTest::LabelledGraph target =
        IsomorphismTest::labelledGraphOf(3, {{0, 1, 1}, {1, 2, 1}});

    VF2 algo(pattern.G, target.G, Semantics::MONOMORPHISM);
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

} // namespace
} // namespace NetworKit
