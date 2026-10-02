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
/// Shared-count updates per worker to aim for, which bounds the overshoot past maxMatches.
constexpr count PublishRounds = 8;

/// Section 5.2.2 of Kimmig, Meyerhenke and Strash finds 4 best.
constexpr count TaskGroupSize = 4;
constexpr count StealAttempts = 4;

/// Thieves use this, so that a steal never waits for the victim.
bool tryLockQueue(std::atomic<bool> &flag) {
    return !flag.exchange(true, std::memory_order_acquire);
}

void lockQueue(std::atomic<bool> &flag) {
    while (!tryLockQueue(flag))
        std::this_thread::yield();
}

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
 * The worker pool of @ref ParallelRI. Every worker walks the search tree depth first on its private
 * `Worker::states` and moves its oldest states to `Worker::stealable`. An idle worker reclaims its
 * own published states first and then steals the shallowest states of another worker. A token
 * passed around the ring of workers detects termination; see @ref passToken().
 */
class ParallelRIImpl {

public:
    /// Passes a match to the user's callback. Returns true if a callback received the match.
    using Deliver = std::function<bool(index, const Match &)>;

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
          publishInterval(
              maxMatches == 0
                  ? 0
                  : std::max<count>(
                        1, std::min<count>(MaxPublishInterval,
                                           maxMatches / (this->numWorkers * PublishRounds)))) {
        for (Worker &worker : workers)
            worker.untilPublish = publishInterval;
    }

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
                // OpenMP may start fewer threads than requested. The ring spans only the started
                // threads, while `workers` is sized for the request.
                activeWorkers.store(
                    std::min<count>(static_cast<count>(omp_get_num_threads()), numWorkers));
                stopOnException([&] { seedRoots(impl); });
            }

            stopOnException([&] { workerLoop(tid, impl); });
        }

        if (failure)
            std::rethrow_exception(failure);
    }

    std::vector<Match> takeMatches() {
        if (!storeMatches)
            return {};

        count total = 0;
        for (const Worker &worker : workers)
            total += worker.buffer.size();

        // The other buffers are appended to the first one.
        std::vector<Match> &merged = workers.front().buffer;
        merged.reserve(total);
        for (index w = 1; w < workers.size(); ++w) {
            std::vector<Match> &buffer = workers[w].buffer;
            merged.insert(merged.end(), std::make_move_iterator(buffer.begin()),
                          std::make_move_iterator(buffer.end()));
            buffer.clear();
            buffer.shrink_to_fit();
        }
        return std::move(merged);
    }

    count matchesFound() const {
        count total = 0;
        for (const Worker &worker : workers)
            total += worker.found;
        return total;
    }

private:
    struct alignas(64) Worker {
        /// Only this worker touches these states.
        std::deque<RIImpl::State> states;
        /// Guarded by `busy`.
        std::deque<RIImpl::State> stealable;
        /// Size of `stealable`, readable without taking `busy`. Exact while `busy` is held.
        std::atomic<count> offered{0};
        std::atomic<bool> busy{false};
        count sinceLastPublish = 0;
        std::vector<Match> buffer;
        count found = 0;
        count untilPublish = 0;
    };

    /// Updates the shared count only every `publishInterval` matches. Returns false once the cap
    /// is reached.
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

    /// An exception must not leave an OpenMP region, so this stores the first one for run() to
    /// rethrow after the join, and stops all workers.
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

    /// Deals the children of the empty mapping round-robin to the workers. This runs inside
    /// `omp single`, so all other workers wait at its barrier.
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

    /// Every state is expanded completely, so `RIImpl::State::nextCandidate` never resumes an
    /// expansion. A position without a parent therefore pushes one state per target node at once.
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
            if (!impl.expand(state, children)) {
                stopped.store(true, std::memory_order_relaxed);
                return;
            }

            for (RIImpl::State &child : children)
                pushLocal(tid, std::move(child));
        }
    }

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

    void pushLocal(index tid, RIImpl::State &&state) {
        Worker &worker = workers[tid];
        worker.states.push_back(std::move(state));
        if (++worker.sinceLastPublish >= TaskGroupSize)
            coalesceIntoTask(tid);
    }

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

    index pickVictim(index thief, count active) const {
        auto &urng = Aux::Random::getURNG();

        // Draw from all workers but one and skip the thief.
        std::uniform_int_distribution<index> pick(0, static_cast<index>(active) - 2);
        const index drawn = pick(urng);
        return drawn < thief ? drawn : drawn + 1;
    }

    /**
     * Passes the termination token on if the idle worker @a tid holds it. The search is over once
     * the token has made a full lap without a successful steal.
     *
     * No lap completes while work remains. A worker with work at the end of a lap was idle at its
     * own hop, so it stole the work afterwards. That steal set `tokenDirty` under the victim's flag
     * before a later hop, which voided the lap. The argument needs sequentially consistent
     * operations on `tokenDirty` and `Worker::offered`.
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

    Aux::SignalHandler *handler;

    Deliver deliver;
    bool storeMatches;
    count maxMatches;
    count numWorkers;

    std::vector<Worker> workers;

    /// Atomic, so that ThreadSanitizer does not report the handover through the OpenMP barrier as
    /// a race.
    std::atomic<count> activeWorkers;

    std::atomic<bool> stopped;
    std::atomic<bool> failed{false};
    std::exception_ptr failure;
    std::atomic<index> tokenHolder;
    /// Set by every successful steal. Voids the current lap of the token.
    std::atomic<bool> tokenDirty;
    std::atomic<count> tokenHops;

    std::atomic<count> published;
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
