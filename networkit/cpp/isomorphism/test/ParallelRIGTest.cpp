/*
 * ParallelRIGTest.cpp
 *
 *  Created on: Aug 24, 2026
 *      Author: Mikhail Kirilin
 */

#include <algorithm>
#include <atomic>
#include <memory>
#include <stdexcept>
#include <vector>

#include <omp.h>
#include <gtest/gtest.h>

#include <networkit/GlobalState.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/isomorphism/ParallelRI.hpp>
#include <networkit/isomorphism/RI.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"

namespace NetworKit {

using IsomorphismTest::graphOf;
using IsomorphismTest::Match;
using IsomorphismTest::sortMatches;
using Semantics = SubgraphIsomorphism::Semantics;

namespace {

/// Builds a ParallelRI of the given variant in the shape the shared harness wants.
template <typename Variant>
auto parallelFactory(Variant variant) {
    return [variant](const Graph &pattern, const Graph &target, Semantics semantics,
                     count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new ParallelRI(pattern, target, variant, semantics, maxMatches));
    };
}

/// The graph the agreement tests run on. Small enough to read in milliseconds, real enough that
/// the search tree is as lopsided as the work stealing is meant to cope with - which the 4-node
/// corpus in the shared harness is far too regular to be.
Graph karate() {
    METISGraphReader reader;
    return reader.read("input/karate.graph");
}

/// A path of @a n nodes. Lengthening it is the cheapest way to buy a deeper and much wider search
/// tree without changing the target: in karate a 5-path occurs 22 064 times and a 7-path 326 328.
Graph path(count n) {
    Graph G(n);
    for (node u = 0; u + 1 < n; ++u)
        G.addEdge(u, u + 1);
    return G;
}

Graph triangle() {
    return graphOf(3, {{0, 1}, {1, 2}, {2, 0}});
}

std::vector<int> workerCounts() {
    std::vector<int> counts;
    for (int workers = 1; workers <= 2 * omp_get_num_procs(); workers *= 2)
        counts.push_back(workers);
    return counts;
}

/// Sorted matches from the sequential search, which is the answer ParallelRI has to reproduce.
std::vector<Match> sequentialMatches(const Graph &pattern, const Graph &target, Semantics semantics,
                                     RI::Variant variant) {
    RI algo(pattern, target, variant, semantics, 0);
    algo.run();
    std::vector<Match> matches = algo.getMatches();
    sortMatches(matches);
    return matches;
}

/// Sorted matches from the parallel search, at whatever the current thread setting is.
std::vector<Match> parallelMatches(const Graph &pattern, const Graph &target, Semantics semantics,
                                   RI::Variant variant) {
    ParallelRI algo(pattern, target, variant, semantics, 0);
    algo.run();
    std::vector<Match> matches = algo.getMatches();
    sortMatches(matches);
    return matches;
}

/// Runs @a algo without storing matches and returns how many it found.
template <typename Algo>
count countOnly(Algo &&algo) {
    algo.setStoreMatches(false);
    algo.run();
    return algo.numberOfMatches();
}

} // namespace

/**
 * Parameterised over the variant, exactly as RIGTest is, so RI-Ds gets the same parallel exercise
 * plain RI does - including the domains, which are built once by the driver and then read by every
 * worker at the same time, a sharing pattern nothing on the sequential path exercises.
 *
 * A hung test here is a termination bug rather than a slow test: the failure mode of the token
 * ring is a worker that never notices the search is over, not a wrong answer.
 */
class ParallelRIGTest : public testing::TestWithParam<RI::Variant> {

protected:
    void SetUp() override { threadsBefore = Aux::getMaxNumberOfThreads(); }

    /// Both of these are process-*global*. A failing ASSERT_* returns from a test early, so
    /// restoring them here rather than at the end of each test is what keeps one failure from
    /// quietly changing how every later test in the binary behaves.
    void TearDown() override {
        GlobalState::setReceivedSIGINT(false);
        Aux::setNumberOfThreads(threadsBefore);
    }

private:
    int threadsBefore = 1;
};

INSTANTIATE_TEST_SUITE_P(Variants, ParallelRIGTest,
                         testing::Values(RI::Variant::RI, RI::Variant::RI_DS));

// -------------------------------------------------------------------------------------------
// The three assertions the shared harness offers
// -------------------------------------------------------------------------------------------

TEST_P(ParallelRIGTest, testAgreesWithTheReference) {

    // The match set, the match cap, and the three callback forms. The last of those is the first
    // real exercise of SubgraphIsomorphism::invokeCallback()'s mutex: a serial MatchCallback must
    // see every match exactly once even though several workers produce them.
    const auto make = parallelFactory(GetParam());

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
    IsomorphismTest::expectCallbackFormsAgree(make);
}

// -------------------------------------------------------------------------------------------
// Agreement with the sequential search, on a graph big enough for the workers to overlap
// -------------------------------------------------------------------------------------------

/**
 * The property that catches lost and duplicated work: the answer must not depend on how many
 * workers produced it.
 *
 * The comparison is element by element rather than by count, because a count alone would let a
 * lost match and a duplicated one cancel each other out - exactly the shape a stealing bug takes.
 *
 * One worker is not a special case in the implementation - it walks the queues, the coalescing and
 * the token ring like any other count - so this really does compare the same machinery at
 * different degrees of contention.
 */
TEST_P(ParallelRIGTest, testAnswerDoesNotDependOnWorkerCount) {
    const Graph target = karate();
    const Graph pattern = path(5);

    const std::vector<Match> expected =
        sequentialMatches(pattern, target, Semantics::MONOMORPHISM, GetParam());
    ASSERT_EQ(expected.size(), 22064u) << "a change here would quietly weaken every case below";

    for (const int workers : workerCounts()) {
        Aux::setNumberOfThreads(workers);

        ParallelRI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
        ASSERT_EQ(algo.numberOfWorkers(), static_cast<count>(workers))
            << "numberOfWorkers() has to follow Aux::setNumberOfThreads()";
        algo.run();

        std::vector<Match> actual = algo.getMatches();
        sortMatches(actual);
        EXPECT_EQ(actual, expected) << "workers: " << workers;
        EXPECT_EQ(algo.numberOfMatches(), expected.size()) << "workers: " << workers;
    }
}

/**
 * The same question again, on the input shape RI-DS-SI-FC exists for: one pattern node whose
 * domain has been narrowed to a single target node.
 *
 * Nothing else in either test file drives a singleton domain through the parallel path, and a
 * singleton is what all three of the RI-Ds improvements turn on - it opens the matching order, it
 * is what forward checking strikes out of every other domain, and it is the reason the order
 * differs from plain RI's. Since the preprocessing result is built once and read concurrently by
 * every worker, a worker reading a domain the driver got wrong would show up here as a wrong
 * answer rather than as a crash.
 *
 * The reference is sequential **plain RI**, not sequential RI-Ds: an independent search that has
 * no domains at all, so agreement cannot come from both sides making the same mistake.
 */
TEST_P(ParallelRIGTest, testSingletonDomainAgreesAtEveryWorkerCount) {
    constexpr node anchor = 0;
    constexpr index rareLabel = 7;

    const Graph target = karate();
    const Graph pattern = path(4);

    // Exactly one target node carries the rare label, and the pattern node that requires it is an
    // endpoint - the position plain RI's degree-first rule would order last.
    std::vector<index> targetNodeLabels(target.upperNodeIdBound(), 0);
    targetNodeLabels[anchor] = rareLabel;
    const std::vector<index> patternNodeLabels{none, none, none, rareLabel};

    RI reference(pattern, target, RI::Variant::RI, Semantics::MONOMORPHISM, 0);
    reference.setNodeLabels(patternNodeLabels, targetNodeLabels);
    reference.run();
    std::vector<Match> expected = reference.getMatches();
    sortMatches(expected);
    ASSERT_FALSE(expected.empty());

    // Every match has to put the labelled pattern node on the one target node that can host it.
    for (const Match &match : expected)
        ASSERT_EQ(match[3], anchor);

    for (const int workers : workerCounts()) {
        Aux::setNumberOfThreads(workers);

        ParallelRI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
        algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
        algo.run();

        std::vector<Match> actual = algo.getMatches();
        sortMatches(actual);
        EXPECT_EQ(actual, expected) << "workers: " << workers;
    }
}

/**
 * Duplication and loss told apart, which comparing sets cannot do, and the worker id contract that
 * makes telling them apart possible at all.
 *
 * The per-worker slots are merged into one sequence and sorted, so the comparison is between
 * multisets: a match reported twice makes the sequence longer, a match lost makes it shorter, and
 * either way the two sequences stop being equal.
 *
 * Sizing those slots from numberOfWorkers() is the whole reason that accessor is public, and the
 * reason the worker count is asked for once inside run() rather than re-read per match. A worker
 * id outside [0, numberOfWorkers()) would make every documented per-worker accumulator write out
 * of bounds, so the id is checked on every match rather than assumed. The thread count is pinned
 * first, so the number the accessor reports is one this test chose.
 */
TEST_P(ParallelRIGTest, testEveryMatchIsReportedExactlyOnce) {
    Aux::setNumberOfThreads(4);

    const Graph target = karate();
    const Graph pattern = triangle();

    const std::vector<Match> expected =
        sequentialMatches(pattern, target, Semantics::MONOMORPHISM, GetParam());
    ASSERT_EQ(expected.size(), 270u) << "karate's triangle count is what pins this case";

    ParallelRI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
    ASSERT_EQ(algo.numberOfWorkers(), 4u);

    std::vector<std::vector<Match>> perWorker(algo.numberOfWorkers());
    algo.setCallback([&](index tid, const Match &match) {
        ASSERT_LT(tid, perWorker.size()) << "a worker id outside [0, numberOfWorkers()) would make "
                                            "every documented per-worker accumulator unsafe";
        perWorker[tid].push_back(match);
    });
    algo.run();

    std::vector<Match> collected;
    for (const std::vector<Match> &slot : perWorker)
        collected.insert(collected.end(), slot.begin(), slot.end());

    sortMatches(collected);
    EXPECT_EQ(collected, expected);
    EXPECT_EQ(algo.numberOfMatches(), expected.size());
}

/**
 * A pattern in two components, so some position has no parent and draws candidates from the whole
 * target rather than from one neighbourhood.
 *
 * That parentless path is also where RI-Ds uses a domain unconditionally, so this is the case in
 * which the two variants do the most different things while having to agree.
 */
TEST_P(ParallelRIGTest, testDisconnectedPatternIsSearchedInParallel) {
    const Graph target = karate();
    const Graph pattern = graphOf(4, {{0, 1}, {2, 3}});

    const std::vector<Match> expected =
        sequentialMatches(pattern, target, Semantics::MONOMORPHISM, GetParam());
    ASSERT_GT(expected.size(), 1000u);

    EXPECT_EQ(parallelMatches(pattern, target, Semantics::MONOMORPHISM, GetParam()), expected);
}

// -------------------------------------------------------------------------------------------
// The inputs where the pool has to stop without ever having had anything to do
// -------------------------------------------------------------------------------------------

TEST_P(ParallelRIGTest, testDegenerateInputs) {
    // Every graph here is a named local, never a temporary: SubgraphIsomorphism holds both graphs
    // by pointer so that a caller can mutate them between runs, which makes `ParallelRI algo(f(),
    // g(), ...)` a dangling read the moment the constructor returns.
    const Graph k4 = graphOf(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
    const Graph nothing(0);
    const Graph oneEdge = graphOf(2, {{0, 1}});
    const Graph path3 = graphOf(3, {{0, 1}, {1, 2}});
    const Graph tri = triangle();

    // An empty pattern has exactly one match - the empty mapping - and it is reported during
    // seeding, before any worker has a queue to drain. No bail-out may swallow it.
    ParallelRI empty(nothing, k4, GetParam(), Semantics::INDUCED, 0);
    empty.run();
    ASSERT_TRUE(empty.hasFinished());
    EXPECT_EQ(empty.numberOfMatches(), 1u);
    ASSERT_EQ(empty.getMatches().size(), 1u);
    EXPECT_TRUE(empty.getMatches().front().empty());

    // More pattern nodes than target nodes: nothing can match, and run() has to say so through
    // patternCannotFit() rather than by starting a pool with an empty seed set.
    ParallelRI tooBig(k4, oneEdge, GetParam(), Semantics::INDUCED, 0);
    tooBig.run();
    ASSERT_TRUE(tooBig.hasFinished());
    EXPECT_EQ(tooBig.numberOfMatches(), 0u);
    EXPECT_FALSE(tooBig.hasMatch());

    // The other way to find nothing: the shape passes every cheap bail-out, so the pool really
    // runs and every branch dies. This is the one that hangs if termination is wrong.
    ParallelRI noMatch(tri, path3, GetParam(), Semantics::MONOMORPHISM, 0);
    noMatch.run();
    ASSERT_TRUE(noMatch.hasFinished());
    EXPECT_EQ(noMatch.numberOfMatches(), 0u);

    // A cap of zero means "no limit", not "no matches" - the same reading RI gives it.
    ParallelRI uncapped(tri, k4, GetParam(), Semantics::MONOMORPHISM, 0);
    uncapped.run();
    EXPECT_EQ(uncapped.numberOfMatches(), 24u);
}

// -------------------------------------------------------------------------------------------
// Interruption and throwing callbacks, which have to unwind every worker and then throw once
// -------------------------------------------------------------------------------------------

/**
 * CTRL+C partway through a long enumeration.
 *
 * The workers may only poll the non-throwing isRunning(), because an exception leaving an OpenMP
 * structured block is undefined behaviour; the InterruptException comes from the single
 * assureRunning() after the join. What this test really proves is that the pool *unwinds* - if any
 * worker kept spinning in the token ring, the region would never join and the binary would hang
 * here rather than fail.
 */
TEST_P(ParallelRIGTest, testInterruptStopsEveryWorker) {
    const Graph target = karate();
    const Graph pattern = path(5);

    ParallelRI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);

    // Several workers hit this at once, so the counter has to be atomic; a user's callback would
    // not need one unless it kept state of its own.
    std::atomic<count> delivered{0};
    algo.setCallback([&](index, const Match &) {
        if (delivered.fetch_add(1) + 1 == 5)
            GlobalState::setReceivedSIGINT(true);
    });

    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);

    EXPECT_FALSE(algo.hasFinished()) << "an interrupted run must not count as finished";
    EXPECT_THROW(algo.numberOfMatches(), std::runtime_error);
    EXPECT_GE(delivered.load(), 5u) << "matches already handed over cannot be taken back";

    // And nothing is poisoned: a clean second search gives the whole answer.
    const std::vector<Match> expected =
        sequentialMatches(pattern, target, Semantics::MONOMORPHISM, GetParam());
    EXPECT_EQ(parallelMatches(pattern, target, Semantics::MONOMORPHISM, GetParam()), expected);
}

/**
 * A callback that throws on every match it is handed, in both callback forms.
 *
 * Several workers can reach the callback at once, and each of them throws. An exception that
 * escaped the OpenMP region would call std::terminate and take the whole test binary down. The
 * first exception has to come out of run() instead, after every worker has unwound. The serial
 * form also has to release its lock on the way out, or the next worker would block on it forever
 * and the region would never join.
 */
TEST_P(ParallelRIGTest, testThrowingCallbackStopsEveryWorker) {
    // Local, so that no other runtime_error out of run() can satisfy EXPECT_THROW.
    struct CallbackFailure : std::runtime_error {
        CallbackFailure() : std::runtime_error("callback failed") {}
    };

    const Graph target = karate();
    const Graph pattern = path(5);
    const count expected = countOnly(RI(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0));

    for (const bool parallelForm : {true, false}) {
        SCOPED_TRACE(parallelForm ? "ParallelMatchCallback" : "MatchCallback");

        ParallelRI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
        std::atomic<count> delivered{0};
        if (parallelForm)
            algo.setCallback([&](index, const Match &) {
                delivered.fetch_add(1);
                throw CallbackFailure();
            });
        else
            algo.setCallback([&](const Match &) {
                delivered.fetch_add(1);
                throw CallbackFailure();
            });

        EXPECT_THROW(algo.run(), CallbackFailure);
        EXPECT_FALSE(algo.hasFinished()) << "a run that threw must not count as finished";

        // Every call throws, so a worker that kept searching after its own exception would call
        // back a second time and push this past one call per worker.
        EXPECT_GE(delivered.load(), 1u);
        EXPECT_LE(delivered.load(), algo.numberOfWorkers());

        // And nothing is poisoned: the same object, given a callback that does not throw, runs to
        // the end and finds everything.
        std::atomic<count> counted{0};
        algo.setCallback([&](index, const Match &) { counted.fetch_add(1); });
        algo.run();
        EXPECT_TRUE(algo.hasFinished());
        EXPECT_EQ(counted.load(), expected);
    }
}

} // namespace NetworKit
