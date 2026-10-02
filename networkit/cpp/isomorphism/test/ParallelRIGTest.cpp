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

template <typename Variant>
auto parallelFactory(Variant variant) {
    return [variant](const Graph &pattern, const Graph &target, Semantics semantics,
                     count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new ParallelRI(pattern, target, variant, semantics, maxMatches));
    };
}

/// Unlike the corpus, karate gives a search tree lopsided enough to exercise the work stealing.
Graph karate() {
    METISGraphReader reader;
    return reader.read("input/karate.graph");
}

/// In karate, a 5-path occurs 22 064 times and a 7-path 326 328 times.
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

std::vector<Match> sequentialMatches(const Graph &pattern, const Graph &target, Semantics semantics,
                                     RI::Variant variant) {
    RI algo(pattern, target, variant, semantics, 0);
    algo.run();
    std::vector<Match> matches = algo.getMatches();
    sortMatches(matches);
    return matches;
}

std::vector<Match> parallelMatches(const Graph &pattern, const Graph &target, Semantics semantics,
                                   RI::Variant variant) {
    ParallelRI algo(pattern, target, variant, semantics, 0);
    algo.run();
    std::vector<Match> matches = algo.getMatches();
    sortMatches(matches);
    return matches;
}

template <typename Algo>
count countOnly(Algo &&algo) {
    algo.setStoreMatches(false);
    algo.run();
    return algo.numberOfMatches();
}

} // namespace

/// A hanging test points to a termination bug: a worker that never notices that the search is over.
class ParallelRIGTest : public testing::TestWithParam<RI::Variant> {

protected:
    void SetUp() override { threadsBefore = Aux::getMaxNumberOfThreads(); }

    /// Restores the process-global state, also after a failing ASSERT_* returned early.
    void TearDown() override {
        GlobalState::setReceivedSIGINT(false);
        Aux::setNumberOfThreads(threadsBefore);
    }

private:
    int threadsBefore = 1;
};

INSTANTIATE_TEST_SUITE_P(Variants, ParallelRIGTest,
                         testing::Values(RI::Variant::RI, RI::Variant::RI_DS));

TEST_P(ParallelRIGTest, testAgreesWithTheReference) {
    const auto make = parallelFactory(GetParam());

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
    IsomorphismTest::expectCallbackFormsAgree(make);
}

/// Compares the matches element by element, since in a count a lost and a duplicated match would
/// cancel out.
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

/// The workers publish their counts in batches and may overshoot the cap, which the stored matches
/// must not show.
TEST_P(ParallelRIGTest, testCapHoldsAtEveryWorkerCount) {
    constexpr count cap = 1000;
    const Graph target = karate();
    const Graph pattern = path(5);

    const std::vector<Match> all =
        sequentialMatches(pattern, target, Semantics::MONOMORPHISM, GetParam());
    ASSERT_GT(all.size(), cap);

    for (const int workers : workerCounts()) {
        Aux::setNumberOfThreads(workers);

        ParallelRI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, cap);
        algo.run();

        EXPECT_EQ(algo.numberOfMatches(), cap) << "workers: " << workers;

        std::vector<Match> actual = algo.getMatches();
        sortMatches(actual);
        EXPECT_EQ(actual.size(), cap) << "workers: " << workers;
        EXPECT_TRUE(std::includes(all.begin(), all.end(), actual.begin(), actual.end()))
            << "workers: " << workers << " - a stored match is invalid or duplicated";
    }
}

/**
 * A singleton domain drives all three RI-DS rules, and all workers read the domains concurrently.
 * The reference is sequential plain RI, which has no domains and so cannot share a domain bug.
 */
TEST_P(ParallelRIGTest, testSingletonDomainAgreesAtEveryWorkerCount) {
    constexpr node anchor = 0;
    constexpr index rareLabel = 7;

    const Graph target = karate();
    const Graph pattern = path(4);

    // Plain RI would order the pattern endpoint that requires the rare label last.
    std::vector<index> targetNodeLabels(target.upperNodeIdBound(), 0);
    targetNodeLabels[anchor] = rareLabel;
    const std::vector<index> patternNodeLabels{none, none, none, rareLabel};

    RI reference(pattern, target, RI::Variant::RI, Semantics::MONOMORPHISM, 0);
    reference.setNodeLabels(patternNodeLabels, targetNodeLabels);
    reference.run();
    std::vector<Match> expected = reference.getMatches();
    sortMatches(expected);
    ASSERT_FALSE(expected.empty());

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

/// The merged per-worker slots form a multiset, so a duplicated or a lost match changes the sorted
/// sequence. Every worker id must lie in [0, numberOfWorkers()).
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

/// A position without a parent draws its candidates from the whole target, or from the domain under
/// RI-DS.
TEST_P(ParallelRIGTest, testDisconnectedPatternIsSearchedInParallel) {
    const Graph target = karate();
    const Graph pattern = graphOf(4, {{0, 1}, {2, 3}});

    const std::vector<Match> expected =
        sequentialMatches(pattern, target, Semantics::MONOMORPHISM, GetParam());
    ASSERT_GT(expected.size(), 1000u);

    EXPECT_EQ(parallelMatches(pattern, target, Semantics::MONOMORPHISM, GetParam()), expected);
}

TEST_P(ParallelRIGTest, testDegenerateInputs) {
    // SubgraphIsomorphism holds the graphs by pointer, so they must not be temporaries.
    const Graph k4 = graphOf(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
    const Graph nothing(0);
    const Graph oneEdge = graphOf(2, {{0, 1}});
    const Graph path3 = graphOf(3, {{0, 1}, {1, 2}});
    const Graph tri = triangle();

    // Seeding reports the only match of an empty pattern.
    ParallelRI empty(nothing, k4, GetParam(), Semantics::INDUCED, 0);
    empty.run();
    ASSERT_TRUE(empty.hasFinished());
    EXPECT_EQ(empty.numberOfMatches(), 1u);
    ASSERT_EQ(empty.getMatches().size(), 1u);
    EXPECT_TRUE(empty.getMatches().front().empty());

    // The empty mapping alone reaches a cap of one, so seeding stops the search.
    ParallelRI emptyCapped(nothing, k4, GetParam(), Semantics::INDUCED, 1);
    emptyCapped.run();
    ASSERT_TRUE(emptyCapped.hasFinished());
    EXPECT_EQ(emptyCapped.numberOfMatches(), 1u);

    ParallelRI tooBig(k4, oneEdge, GetParam(), Semantics::INDUCED, 0);
    tooBig.run();
    ASSERT_TRUE(tooBig.hasFinished());
    EXPECT_EQ(tooBig.numberOfMatches(), 0u);
    EXPECT_FALSE(tooBig.hasMatch());

    // The pool runs and every branch dies, which hangs if termination is wrong.
    ParallelRI noMatch(tri, path3, GetParam(), Semantics::MONOMORPHISM, 0);
    noMatch.run();
    ASSERT_TRUE(noMatch.hasFinished());
    EXPECT_EQ(noMatch.numberOfMatches(), 0u);

    ParallelRI uncapped(tri, k4, GetParam(), Semantics::MONOMORPHISM, 0);
    uncapped.run();
    EXPECT_EQ(uncapped.numberOfMatches(), 24u);
}

/// The binary hangs here instead of failing if a worker keeps spinning after CTRL+C.
TEST_P(ParallelRIGTest, testInterruptStopsEveryWorker) {
    const Graph target = karate();
    const Graph pattern = path(5);

    ParallelRI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);

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

    const std::vector<Match> expected =
        sequentialMatches(pattern, target, Semantics::MONOMORPHISM, GetParam());
    EXPECT_EQ(parallelMatches(pattern, target, Semantics::MONOMORPHISM, GetParam()), expected);
}

/// An exception that escaped the OpenMP region would terminate the binary. The serial callback form
/// must also release its lock, or the next worker would block forever.
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

        // A worker that kept searching after its own exception would call back a second time.
        EXPECT_GE(delivered.load(), 1u);
        EXPECT_LE(delivered.load(), algo.numberOfWorkers());

        std::atomic<count> counted{0};
        algo.setCallback([&](index, const Match &) { counted.fetch_add(1); });
        algo.run();
        EXPECT_TRUE(algo.hasFinished());
        EXPECT_EQ(counted.load(), expected);
    }
}

} // namespace NetworKit
