#include <algorithm>
#include <atomic>
#include <deque>
#include <exception>
#include <functional>
#include <iterator>
#include <random>
#include <thread>
#include <utility>
#include <vector>

#include <omp.h>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/isomorphism/ParallelRI.hpp>

#include "RIImpl.hpp"
#include "SearchGraph.hpp"

namespace NetworKit {

namespace {

using Match = SubgraphIsomorphism::Match;

using IsomorphismDetails::RIImpl;
using IsomorphismDetails::SearchGraph;

/// Maximum number of matches a worker records between two updates of the shared count.
constexpr count MaxPublishInterval = 64;
/// Number of shared-count updates per worker to aim for, which bounds the overshoot past
/// maxMatches.
constexpr count PublishRounds = 8;

/// Number of states a worker publishes for stealing at once. Section V-B2 of the paper finds 4
/// best.
constexpr count TaskGroupSize = 4;
/// Number of victims a thief tries before it joins the termination protocol.
constexpr count StealAttempts = 4;

/// Tries to take a worker's queue flag without waiting. Thieves use this, so that a steal never
/// waits for the victim.
bool tryLockQueue(std::atomic<bool> &flag) {
    return !flag.exchange(true, std::memory_order_acquire);
}

/// Takes a worker's queue flag, spinning until it is free. The owner uses this.
void lockQueue(std::atomic<bool> &flag) {
    while (!tryLockQueue(flag))
        std::this_thread::yield();
}

/// Releases a queue flag taken by @ref lockQueue() or @ref tryLockQueue() when the scope is left.
class QueueUnlock {

public:
    explicit QueueUnlock(std::atomic<bool> &flag) : flag(&flag) {}
    QueueUnlock(const QueueUnlock &) = delete;
    QueueUnlock &operator=(const QueueUnlock &) = delete;
    ~QueueUnlock() { flag->store(false, std::memory_order_release); }

private:
    std::atomic<bool> *flag;
};

/**
 * The worker pool of @ref ParallelRI. Every worker owns an RIImpl and a private queue of states.
 *
 * A worker pushes and pops states at the back of `Worker::states`, which no other thread touches,
 * so it walks the search tree depth first without synchronization. When enough states have
 * accumulated, it moves the oldest ones to `Worker::stealable`, which `Worker::busy` guards. A
 * worker without states first reclaims its own `stealable` states and then steals from the front
 * of another worker's `stealable`, where the shallowest states with the largest subtrees are. A
 * token passed around the ring of workers detects termination; see @ref passToken().
 *
 * Every worker counts and stores its matches in its own slot, and ParallelRI::run() merges them
 * after the join. Only the user's callback is shared, through @a deliver.
 */
class ParallelRIImpl {

public:
    /// Passes a match to the user's callback. Thread-safe. Returns true if a callback received the
    /// match, in which case the match must not be stored.
    using Deliver = std::function<bool(index, const Match &)>;

    /**
     * @param patternGraph Snapshot of the pattern, shared read-only.
     * @param targetGraph Snapshot of the target, shared read-only.
     * @param patternNodeLabels Empty if the search is unlabelled.
     * @param targetNodeLabels Empty if the search is unlabelled.
     * @param ordering Matching order, shared read-only.
     * @param domains RI-DS domains, shared read-only. Empty under plain RI.
     * @param semantics Whether matches must be induced.
     * @param handler Shared by all workers. Only isRunning() may be called on it inside the
     * parallel region.
     * @param deliver Passes a match to the user's callback. Called concurrently.
     * @param storeMatches Whether matches have to be stored rather than only counted.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     * @param numWorkers Number of workers, as reported by SubgraphIsomorphism::numberOfWorkers().
     */
    ParallelRIImpl(const SearchGraph &patternGraph, const SearchGraph &targetGraph,
                   const std::vector<index> &patternNodeLabels,
                   const std::vector<index> &targetNodeLabels, const RIImpl::Ordering &ordering,
                   const RIImpl::Domains &domains, SubgraphIsomorphism::Semantics semantics,
                   Aux::SignalHandler &handler, Deliver deliver, bool storeMatches,
                   count maxMatches, count numWorkers)
        : patternGraph(&patternGraph), targetGraph(&targetGraph),
          patternNodeLabels(&patternNodeLabels), targetNodeLabels(&targetNodeLabels),
          ordering(&ordering), domains(&domains), semantics(semantics), handler(&handler),
          deliver(std::move(deliver)), storeMatches(storeMatches), maxMatches(maxMatches),
          numWorkers(numWorkers == 0 ? 1 : numWorkers),
          // Worker is neither copyable nor movable, so the vector is sized on construction.
          workers(this->numWorkers), activeWorkers(this->numWorkers), stopped(false),
          tokenHolder(0), tokenDirty(false), tokenHops(0), published(0),
          // Publish the local count rarely, but often enough to overshoot maxMatches only slightly.
          publishInterval(
              maxMatches == 0
                  ? 0
                  : std::max<count>(
                        1, std::min<count>(MaxPublishInterval,
                                           maxMatches / (this->numWorkers * PublishRounds)))) {
        for (Worker &worker : workers)
            worker.untilPublish = publishInterval;
    }

    /**
     * Runs the parallel search. If a worker caught an exception, rethrows the first one after all
     * workers have joined.
     */
    void run() {
        // The early exits of RIImpl::run(). An empty pattern has one match, the empty mapping.
        if (!ordering->order.empty()
            && (RIImpl::patternCannotFit(*patternGraph, *targetGraph) || domains->anyEmpty))
            return;

#pragma omp parallel num_threads(static_cast<int>(numWorkers))
        {
            const index tid = static_cast<index>(omp_get_thread_num());

            RIImpl impl(*patternGraph, *targetGraph, *patternNodeLabels, *targetNodeLabels,
                        *ordering, *domains, semantics, *handler,
                        [this, tid](const Match &match) { return recordMatch(tid, match); });

#pragma omp single
            {
                // OpenMP may start fewer threads than requested, but never more. The ring spans
                // only the started threads, while `workers` is sized for the request.
                activeWorkers.store(
                    std::min<count>(static_cast<count>(omp_get_num_threads()), numWorkers));
                stopOnException([&] { seedRoots(impl); });
            }
            // The barrier at the end of the single region lets the seeding finish before any
            // worker starts.

            stopOnException([&] { workerLoop(tid, impl); });
        }

        if (failure)
            std::rethrow_exception(failure);
    }

    /**
     * @return the concatenated matches of all workers. Empty if no matches were stored. Call this
     * once, after @ref run().
     */
    std::vector<Match> takeMatches() {
        std::vector<Match> merged;
        if (!storeMatches)
            return merged;

        count total = 0;
        for (const Worker &worker : workers)
            total += worker.buffer.size();

        merged.reserve(total);
        for (Worker &worker : workers) {
            merged.insert(merged.end(), std::make_move_iterator(worker.buffer.begin()),
                          std::make_move_iterator(worker.buffer.end()));
            worker.buffer.clear();
            worker.buffer.shrink_to_fit();
        }
        return merged;
    }

    /// @return the total number of matches found, stored or not. Call this after @ref run().
    count matchesFound() const {
        count total = 0;
        for (const Worker &worker : workers)
            total += worker.found;
        return total;
    }

private:
    /// The state of one worker, padded to its own cache line.
    struct alignas(64) Worker {
        /// States that only this worker touches, pushed and popped at the back.
        std::deque<RIImpl::State> states;
        /// States published for stealing. Guarded by `busy`.
        std::deque<RIImpl::State> stealable;
        /// Size of `stealable`, readable without taking `busy`. Exact while `busy` is held.
        std::atomic<count> offered{0};
        /// Guards `stealable`.
        std::atomic<bool> busy{false};
        /// Number of states pushed since the last publication.
        count sinceLastPublish = 0;
        /// Matches found and stored by this worker.
        std::vector<Match> buffer;
        /// Number of matches found by this worker, stored or not.
        count found = 0;
        /// Matches left until the next update of `published`. Only used if maxMatches != 0.
        count untilPublish = 0;
    };

    /**
     * Records a match found by worker @a tid. With a cap on the number of matches, the shared
     * count is only updated every `publishInterval` matches.
     *
     * @return false once the cap is reached.
     */
    bool recordMatch(index tid, const Match &match) {
        Worker &worker = workers[tid];
        ++worker.found;

        if (!deliver(tid, match) && storeMatches)
            worker.buffer.push_back(match);

        if (maxMatches == 0)
            return true;

        // A single worker applies the cap exactly.
        if (activeWorkers.load(std::memory_order_relaxed) == 1)
            return worker.found < maxMatches;

        if (--worker.untilPublish == 0) {
            worker.untilPublish = publishInterval;
            const count total =
                published.fetch_add(publishInterval, std::memory_order_relaxed) + publishInterval;
            if (total >= maxMatches) {
                stopped.store(true, std::memory_order_relaxed);
                return false;
            }
        }

        return !stopped.load(std::memory_order_relaxed);
    }

    /**
     * Runs @a step. If it throws, stores the first exception, which run() rethrows after the join,
     * and stops all workers. An exception must not leave an OpenMP region.
     */
    template <typename Step>
    void stopOnException(Step &&step) {
        try {
            step();
        } catch (...) {
            if (!failed.exchange(true))
                failure = std::current_exception();
            stopped.store(true, std::memory_order_relaxed);
        }
    }

    /**
     * Expands the empty mapping and deals the resulting states round-robin to the workers. For an
     * empty pattern, this reports the only match. Called inside `omp single`, while all other
     * workers wait at its barrier.
     *
     * @param impl The RIImpl of the seeding thread.
     */
    void seedRoots(RIImpl &impl) {
        RIImpl::State root = impl.rootState();

        std::vector<RIImpl::State> roots;
        if (!impl.expand(root, roots)) {
            stopped.store(true, std::memory_order_relaxed);
            return;
        }

        const count active = activeWorkers.load();
        index next = 0;
        for (RIImpl::State &state : roots) {
            workers[next].states.push_back(std::move(state));
            next = static_cast<index>((next + 1) % active);
        }
    }

    /**
     * Main loop of worker @a tid: takes a state from its own queue or steals one, expands it and
     * pushes the children. Without a state, it passes the termination token and returns once the
     * search is over.
     *
     * Every state is expanded completely, so `RIImpl::State::nextCandidate` is never used to
     * resume. A position without a parent therefore pushes one state per target node at once.
     */
    void workerLoop(index tid, RIImpl &impl) {
        RIImpl::State state;
        std::vector<RIImpl::State> children;

        while (!stopped.load(std::memory_order_relaxed)) {
            if (!popLocal(tid, state) && !trySteal(tid, state)) {
                passToken(tid);
                if (quiescent())
                    return;
                std::this_thread::yield();
                continue;
            }

            // Only the non-throwing isRunning() may be called inside the parallel region.
            if (!handler->isRunning()) {
                stopped.store(true, std::memory_order_relaxed);
                return;
            }

            children.clear();
            // A state at full depth is reported by expand().
            if (!impl.expand(state, children)) {
                // The cap is reached.
                stopped.store(true, std::memory_order_relaxed);
                return;
            }

            for (RIImpl::State &child : children)
                pushLocal(tid, std::move(child));
        }
    }

    /**
     * Takes the newest state of worker @a tid, reclaiming its published states if its private
     * queue is empty.
     *
     * @return false if the worker has no states left, that is, it is idle.
     */
    bool popLocal(index tid, RIImpl::State &out) {
        Worker &worker = workers[tid];

        if (!worker.states.empty()) {
            out = std::move(worker.states.back());
            worker.states.pop_back();
            return true;
        }

        // Reclaim published states, so that the worker never reports idle while it has work. Only
        // the owner adds to `stealable`, so a zero read here is exact.
        if (worker.offered.load() == 0)
            return false;

        lockQueue(worker.busy);
        const QueueUnlock guard(worker.busy);
        if (worker.stealable.empty())
            return false;

        // Take the newest state and leave the oldest ones to thieves.
        out = std::move(worker.stealable.back());
        worker.stealable.pop_back();
        worker.offered.store(worker.stealable.size());
        return true;
    }

    /// Pushes a state onto the private queue of worker @a tid, and tries to publish a batch every
    /// TaskGroupSize pushes.
    void pushLocal(index tid, RIImpl::State &&state) {
        Worker &worker = workers[tid];
        worker.states.push_back(std::move(state));
        if (++worker.sinceLastPublish >= TaskGroupSize)
            coalesceIntoTask(tid);
    }

    /**
     * Publishes a batch of the oldest states of worker @a tid for stealing, unless states are
     * still on offer or the worker has at most TaskGroupSize states. The worker keeps at least
     * TaskGroupSize states.
     */
    void coalesceIntoTask(index tid) {
        Worker &worker = workers[tid];
        worker.sinceLastPublish = 0;

        if (worker.offered.load() != 0)
            return;

        if (worker.states.size() <= TaskGroupSize)
            return;

        const count batch = std::min<count>(TaskGroupSize, worker.states.size() - TaskGroupSize);

        lockQueue(worker.busy);
        const QueueUnlock guard(worker.busy);
        for (count i = 0; i < batch; ++i) {
            worker.stealable.push_back(std::move(worker.states.front()));
            worker.states.pop_front();
        }
        worker.offered.store(worker.stealable.size());
    }

    /**
     * Tries to steal the oldest published state of another worker.
     *
     * @return false if no state was stolen.
     */
    bool trySteal(index thief, RIImpl::State &out) {
        // pickVictim() needs at least two workers.
        const count active = activeWorkers.load();
        if (active < 2)
            return false;

        for (count attempt = 0; attempt < StealAttempts; ++attempt) {
            Worker &victim = workers[pickVictim(thief, active)];

            if (victim.offered.load() == 0)
                continue;
            if (!tryLockQueue(victim.busy))
                continue;

            const QueueUnlock guard(victim.busy);
            if (victim.stealable.empty())
                continue;

            // Set while holding the victim's flag, which passToken() relies on.
            tokenDirty.store(true);

            out = std::move(victim.stealable.front());
            victim.stealable.pop_front();
            victim.offered.store(victim.stealable.size());
            return true;
        }

        return false;
    }

    /**
     * @param thief The stealing worker, which is never returned.
     * @param active The number of workers; at least 2.
     * @return a uniformly random worker other than @a thief.
     */
    index pickVictim(index thief, count active) const {
        auto &urng = Aux::Random::getURNG();

        // Draw from all workers but one and skip the thief.
        std::uniform_int_distribution<index> pick(0, static_cast<index>(active) - 2);
        const index drawn = pick(urng);
        return drawn < thief ? drawn : drawn + 1;
    }

    /**
     * Passes the termination token on if the idle worker @a tid holds it. The search is over once
     * the token has made a full lap of `activeWorkers` hops without a successful steal.
     *
     * No lap completes while work remains. A worker calls this only after @ref popLocal() found
     * both its queues empty, and every successful steal sets `tokenDirty` while holding the
     * victim's flag. A worker that has work at the end of a lap was idle at its own hop, so it
     * stole the work afterwards. That steal set `tokenDirty` before a later hop of the lap, which
     * then voided the lap. `tokenDirty` and `Worker::offered` use sequentially consistent
     * operations for this argument.
     */
    void passToken(index tid) {
        if (tokenHolder.load() != tid)
            return;

        if (tokenDirty.exchange(false))
            tokenHops.store(0);
        else
            tokenHops.fetch_add(1);

        const count active = activeWorkers.load();
        tokenHolder.store(static_cast<index>((tid + 1) % active));
    }

    /// @return true if the token has made a full lap or the search was stopped.
    bool quiescent() const {
        return stopped.load(std::memory_order_relaxed) || tokenHops.load() >= activeWorkers.load();
    }

    const SearchGraph *patternGraph;
    const SearchGraph *targetGraph;

    const std::vector<index> *patternNodeLabels;
    const std::vector<index> *targetNodeLabels;

    const RIImpl::Ordering *ordering;
    const RIImpl::Domains *domains;

    SubgraphIsomorphism::Semantics semantics;

    /// Only isRunning() may be called on it inside the parallel region.
    Aux::SignalHandler *handler;

    Deliver deliver;
    bool storeMatches;
    /// 0 means no limit.
    count maxMatches;
    count numWorkers;

    std::vector<Worker> workers;

    /// Number of threads that OpenMP actually started, at most `numWorkers`. Atomic, so that
    /// ThreadSanitizer does not report the handover through the OpenMP barrier as a race.
    std::atomic<count> activeWorkers;

    /// Stops all workers once the cap is reached, the search is interrupted or a worker threw.
    std::atomic<bool> stopped;
    /// Set by the first worker that catches an exception.
    std::atomic<bool> failed{false};
    /// The first exception a worker caught. Rethrown by run() after the join.
    std::exception_ptr failure;
    /// The worker that holds the termination token.
    std::atomic<index> tokenHolder;
    /// Set by every successful steal. Voids the current lap of the token.
    std::atomic<bool> tokenDirty;
    /// Hops since the last void. Reaching `activeWorkers` ends the search.
    std::atomic<count> tokenHops;

    /// Matches published by all workers. Only used if maxMatches != 0 and several workers run.
    std::atomic<count> published;
    /// Number of matches a worker records between two updates of `published`.
    count publishInterval;
};

} // namespace

ParallelRI::ParallelRI(const Graph &pattern, const Graph &target, RI::Variant variant,
                       Semantics semantics, count maxMatches)
    : SubgraphIsomorphism(pattern, target, semantics, maxMatches), variant(variant) {}

count ParallelRI::numberOfWorkers() const {
    return static_cast<count>(Aux::getMaxNumberOfThreads());
}

void ParallelRI::run() {
    Aux::SignalHandler handler;

    prepareRun();

    // Read once, so that the search uses the number of workers that numberOfWorkers() reports.
    const count numWorkers = numberOfWorkers();

    // All workers share this setup read-only. prepareRISearch() throws before any worker starts.
    const IsomorphismDetails::RISearchSetup setup = IsomorphismDetails::prepareRISearch(
        *pattern, *target, patternNodeLabels, targetNodeLabels, patternEdgeLabels, targetEdgeLabels,
        variant, "ParallelRI");

    ParallelRIImpl impl(
        setup.patternGraph, setup.targetGraph, patternNodeLabels, targetNodeLabels, setup.ordering,
        setup.domains, semantics, handler,
        [this](index tid, const Match &match) { return invokeCallback(tid, match); },
        storesMatches(), maxMatches, numWorkers);
    impl.run();

    // Throw only after all workers have joined.
    handler.assureRunning();

    finishRun(impl.takeMatches(), impl.matchesFound());
}

} // namespace NetworKit
