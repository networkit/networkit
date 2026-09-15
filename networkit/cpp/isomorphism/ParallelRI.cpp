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

/// Upper bound on how many matches a worker records between two publications of its local count.
constexpr count MaxPublishInterval = 64;
/// How many publication rounds per worker to aim for, bounding the overshoot past maxMatches.
constexpr count PublishRounds = 8;

/// How many states become visible to thieves at a time. Four is what the paper measures as the
/// best trade between the cost of publishing and how long a stolen task keeps a thief busy; its
/// Section V-B2 finds 4 best on both PDBSv1 and GRAEMLIN32 and 16 clearly worse.
constexpr count TaskGroupSize = 4;
/// How many victims a thief tries before giving up and joining the termination protocol.
constexpr count StealAttempts = 4;

/// Take a worker's published-queue flag, or report that somebody else holds it.
///
/// The critical sections guarded by these are a couple of moves long, so spinning beats anything
/// the operating system could offer. The asymmetry between the two entry points is the point: an
/// owner waits for its own flag, a thief never waits at all, so a steal cannot stall behind the
/// victim's expansion.
bool tryLockQueue(std::atomic<bool> &flag) {
    return !flag.exchange(true, std::memory_order_acquire);
}

void lockQueue(std::atomic<bool> &flag) {
    while (!tryLockQueue(flag))
        std::this_thread::yield();
}

/// Releases what @ref lockQueue() or @ref tryLockQueue() took, however the scope is left.
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
 * The worker pool that drives one RIImpl per thread.
 *
 * ## What each worker owns
 *
 * A private queue of partial mappings, a private RIImpl to expand them with, and a private buffer
 * and counter for the matches it finds. None of those are shared, so the common path - take a
 * state off my own queue, expand it, push the children back, record a match - touches no memory
 * any other thread writes to.
 *
 * The only shared things are: the graph snapshots and the matching order, which are read-only
 * after construction; the victims' *published* queues, touched only when stealing; and a handful
 * of atomics used for stopping. That is the whole synchronization surface.
 *
 * ## Where the results go
 *
 * Reporting a match is the one place a worker has something to say to the outside world, and the
 * naive design - one shared result vector behind a lock - would put that lock in the middle of the
 * search. So the accumulation lives here instead, per worker, and ParallelRI::run() merges the
 * buffers once after the join. This class is the only one in the module that needs any of it:
 * VF2, VF3 and RI are sequential and just call SubgraphIsomorphism::reportMatch().
 *
 * The one part that cannot be private is the user's callback, which the caller passes in as
 * @a deliver. It hands a match to whichever callback form was set and handles the serial form's
 * "never entered twice at once" promise itself, so a worker just calls it and moves on.
 *
 * ## Two queues, not one deque with two ends
 *
 * The obvious shape is a single deque with the owner at one end and thieves at the other, and it
 * cannot be built on std::deque: `push_back` and `pop_front` both mutate the container's own
 * bookkeeping, and one of them may reallocate the map array while the other is reading it.
 * Concurrently modifying distinct *elements* of a deque is safe; concurrently modifying the deque
 * is not. So "private" is made literally true by splitting the queue in two.
 *
 * - `Worker::states` is touched by nobody but its owner, ever, and needs no atomic of any kind.
 *   The owner pushes and pops at the **back**, so the walk is depth-first: a child is expanded
 *   before its siblings, and the queue therefore holds at most one level's children per level
 *   below the root, rather than the whole breadth-first frontier. That is *not* the same as
 *   "proportional to the pattern's depth", and the difference matters on one shape in particular:
 *   @ref RIImpl::expand() produces **every** child of a state in one call, so a position with no
 *   parent - the first position, and any position that starts a new component of a disconnected
 *   pattern - hands over one State per existing target node at once. On a large target that is a
 *   queue of `n` full-width mappings. Bounding it would mean taking a state's children in chunks,
 *   which `RIImpl::State::nextCandidate` already allows for and which @ref workerLoop() does not
 *   use today - see the note there.
 * - `Worker::stealable` holds the batches the owner has published, and is the only thing the
 *   `Worker::busy` flag protects. A thief takes from its **front** - the oldest and therefore
 *   shallowest state on offer, which is the one with the most search tree hanging off it - so one
 *   steal buys a lot of work and steals stay rare.
 *
 * The one rule termination depends on: @ref popLocal() reclaims from `stealable` before it reports
 * failure. A worker that still had published states but called itself idle would let the token
 * complete a lap with work outstanding, and the search would end early with a wrong answer.
 */
class ParallelRIImpl {

public:
    /// Hands one match to the user's callback. Thread-safe; returns true if a callback took it,
    /// in which case this class must not store it.
    using Deliver = std::function<bool(index, const Match &)>;

    /**
     * @param patternGraph Shared read-only snapshot of the pattern.
     * @param targetGraph Shared read-only snapshot of the target.
     * @param patternNodeLabels Empty when the search is unlabelled.
     * @param targetNodeLabels Empty when the search is unlabelled.
     * @param ordering Shared read-only matching order, computed once by the caller.
     * @param domains Shared read-only RI-Ds domains, computed once by the caller before the order.
     *        Empty under plain RI, which is how the variant reaches the workers.
     * @param semantics Whether matches must be induced.
     * @param handler Shared by every worker. Only its non-throwing isRunning() may be called
     *        from inside the parallel region.
     * @param deliver Where the user's callback is invoked. Called from several threads at once.
     * @param storeMatches Whether the matches have to be kept, rather than only counted.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     * @param numWorkers How many workers to run, which is what the caller told the user to
     *        expect from SubgraphIsomorphism::numberOfWorkers().
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
          // Worker holds an atomic, so it is neither copyable nor movable and vector::resize()
          // would not compile. vector(size_type) only needs default construction, hence here.
          workers(this->numWorkers), activeWorkers(this->numWorkers), stopped(false),
          tokenHolder(0), tokenDirty(false), tokenHops(0), published(0),
          // Publish often enough that the workers overshoot maxMatches only slightly, but rarely
          // enough that a large enumeration performs no atomic operations worth speaking of.
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
     * Run the parallel search.
     *
     * Everything shared is read-only and already built, so the region below is short: each thread
     * gives itself an RIImpl bound to its own slot, one thread deals out the roots, and then every
     * thread runs the same loop. There is no shortcut for a team of one - a single worker walks
     * the queues, the coalescing and the token ring exactly as sixteen do, which is what makes
     * comparing worker counts a real test and what the paper's own speedup baseline measures.
     *
     * If a worker caught an exception, run() rethrows the first one after every worker has joined.
     * In practice that exception comes from the user's callback; see @ref stopOnException().
     */
    void run() {
        // The bail-outs RIImpl::run() does and the expand() path does not: a pattern with more
        // nodes or a higher maximum degree than the target can match nothing, and under RI-Ds
        // neither can one whose domains left a pattern node with nowhere to go. Guarded on a
        // non-empty order, because an empty pattern has exactly one match - the empty mapping.
        if (!ordering->order.empty()
            && (RIImpl::patternCannotFit(*patternGraph, *targetGraph) || domains->anyEmpty))
            return;

#pragma omp parallel num_threads(static_cast<int>(numWorkers))
        {
            const index tid = static_cast<index>(omp_get_thread_num());

            // One per worker, and cheap: it sizes a buffer and takes two pointers. The domains it
            // reads were built once by the caller, so sixteen workers share one copy of what used
            // to be sixteen scans of the whole target.
            RIImpl impl(*patternGraph, *targetGraph, *patternNodeLabels, *targetNodeLabels,
                        *ordering, *domains, semantics, *handler,
                        [this, tid](const Match &match) { return recordMatch(tid, match); });

#pragma omp single
            {
                // OpenMP may give a smaller team than asked for, and never a larger one - the
                // num_threads clause is an upper bound, which is the whole reason `workers` can be
                // sized once, up front. Everything that has to agree on how many workers there
                // really are reads this, never numWorkers: seeding a queue nobody drains would
                // hang the token ring on a worker that does not exist.
                //
                // The min() clamps the *ring*, and only that. It does not protect the
                // `workers[tid]` indexing in recordMatch(), popLocal() and pushLocal(), which uses
                // the raw thread number: a runtime that broke the "never larger" promise would run
                // past the end of the vector there whatever this line said. Making that safe would
                // mean sizing `workers` from the team, which cannot be done before the region
                // opens. It is the OpenMP contract that holds it up, not this clamp.
                activeWorkers.store(
                    std::min<count>(static_cast<count>(omp_get_num_threads()), numWorkers));
                stopOnException([&] { seedRoots(impl); });
            }
            // The implicit barrier ending the single is what makes writing into other workers'
            // private queues above legal: nobody else has started yet.

            stopOnException([&] { workerLoop(tid, impl); });
        }

        // Every worker has joined, so the exception may finally leave.
        if (failure)
            std::rethrow_exception(failure);
    }

    /**
     * The workers' buffers concatenated. Empty when nothing was stored.
     *
     * Call once, after @ref run() has returned and every worker has joined.
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

    /// How many matches were reported in total, stored or not. Call after @ref run().
    count matchesFound() const {
        count total = 0;
        for (const Worker &worker : workers)
            total += worker.found;
        return total;
    }

private:
    /// One worker's private state. Padded so two workers never share a cache line.
    struct alignas(64) Worker {
        /// Partial mappings only this worker may touch. Pushed and popped at the back, so the walk
        /// is depth-first: at most one expansion's worth of children per level, not the whole
        /// frontier. One expansion can still be large - see the class documentation.
        std::deque<RIImpl::State> states;
        /// The oldest states, published for stealing. Everything here is behind `busy`.
        std::deque<RIImpl::State> stealable;
        /// Size of `stealable`, readable without the lock. This is the paper's work_available: it
        /// lets a thief skip an empty victim, and the owner skip a pointless publication. Exact
        /// under the lock, a hint outside it - every read outside is rechecked once the lock is
        /// held. Sequentially consistent because the termination argument leans on it; see
        /// @ref passToken().
        std::atomic<count> offered{0};
        /// Guards `stealable`. A thief that loses the race moves on rather than waiting.
        std::atomic<bool> busy{false};
        /// How many states have been pushed since the last batch was published for stealing.
        count sinceLastPublish = 0;
        /// Matches this worker found and kept. Merged by takeMatches() after the join.
        std::vector<Match> buffer;
        /// How many matches this worker reported, whether or not they were kept.
        count found = 0;
        /// Counts down to the next publication of `found`; only used when maxMatches != 0.
        count untilPublish = 0;
    };

    /**
     * Record one match found by worker @a tid.
     *
     * This is what each worker's RIImpl reports through. Everything it touches is either that
     * worker's own slot or, when a cap was requested, one atomic that is written only every
     * `publishInterval` matches - so the common case performs no synchronization at all.
     *
     * @return true to keep searching, false once the cap has been reached.
     */
    bool recordMatch(index tid, const Match &match) {
        Worker &worker = workers[tid];
        ++worker.found;

        if (!deliver(tid, match) && storeMatches)
            worker.buffer.push_back(match);

        if (maxMatches == 0)
            return true;

        // A lone worker knows the global count exactly, so it needs no coordination at all. This
        // asks activeWorkers rather than numWorkers, so that a team OpenMP shrank to one still
        // gets the exact cap instead of the publishInterval approximation.
        if (activeWorkers.load(std::memory_order_relaxed) == 1)
            return worker.found < maxMatches;

        // Several workers: tell the others how far we have got, but only every so often, and
        // otherwise just read the flag they set when the cap is reached. The read is a plain load
        // from a cache line nobody writes to, so it is effectively free.
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
     * Run @a step, and turn anything it throws into a stop of the whole search.
     *
     * An exception must not leave an OpenMP structured block. The runtime answers with
     * std::terminate, so a callback that throws would kill the process instead of reaching the
     * caller. A Python callback that raises is the common case: its wrapper rethrows the error as
     * a std::runtime_error, and several workers may do so at once.
     *
     * This function keeps the first exception and sets `stopped`, which makes every other worker
     * unwind; run() rethrows the exception after the join. Any later exception is dropped, because
     * it can only come from a worker that had not seen `stopped` yet. Unwinding out of the middle
     * of an expansion leaves nothing held: the serial callback's mutex sits in a lock_guard, and
     * every queue flag here sits in a QueueUnlock.
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
     * Create the initial states and spread them over the workers.
     *
     * One expansion of the empty mapping gives the ways of placing the first pattern node in the
     * order, which is exactly the paper's initial task set. For an empty pattern this is the whole
     * search, and the match reported here is the empty mapping itself - which is also why the
     * cap has to be honoured on this path.
     *
     * Called from inside `omp single`, so writing into other workers' private queues is legal:
     * they are all waiting at the barrier that ends it.
     *
     * @param impl The seeding thread's own RIImpl, whose reporter is bound to that thread's slot.
     */
    void seedRoots(RIImpl &impl) {
        RIImpl::State root = impl.rootState();

        std::vector<RIImpl::State> roots;
        if (!impl.expand(root, roots)) {
            stopped.store(true, std::memory_order_relaxed);
            return;
        }

        // Round-robin rather than in blocks, so every worker has something to do immediately
        // instead of all of them starting by stealing from worker 0.
        const count active = activeWorkers.load();
        index next = 0;
        for (RIImpl::State &state : roots) {
            workers[next].states.push_back(std::move(state));
            next = static_cast<index>((next + 1) % active);
        }
    }

    /**
     * What one worker does for the whole search.
     *
     * Take a state - my own if I have one, somebody else's if I have not - expand it by one level,
     * push the children back, repeat. With nothing left to take, say so by moving the termination
     * token on, and stop once it has been all the way round with nobody having found work.
     *
     * Note what this does *not* do with @ref RIImpl::State::nextCandidate. Every state reaching
     * @ref RIImpl::expand() here arrives with that field at 0 and is dropped as soon as expand()
     * returns, so the resume point is written and never read: this loop always takes a state's
     * whole child list in one call. That is what makes an expansion of a parentless position hand
     * over one State per target node at once, which the class documentation calls out as the
     * queue's worst case. Feeding a state back onto the queue with its resume point advanced,
     * instead of draining it, is the change that would bound the queue - `expand()` already
     * supports it, and it is the reason the field exists.
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

            // isRunning() and never assureRunning(): an exception leaving an OpenMP structured
            // block is undefined behaviour. ParallelRI::run() throws once, after the join.
            if (!handler->isRunning()) {
                stopped.store(true, std::memory_order_relaxed);
                return;
            }

            children.clear();
            // A state that is already at full depth is reported by expand() itself, so a match
            // costs one queue round trip and every matching rule stays inside RIImpl - which is
            // the module's rule against the two searches drifting apart.
            if (!impl.expand(state, children)) {
                // The cap was reached inside the reporter. Everyone else finds out through
                // `stopped`, which is also what ends their token laps.
                stopped.store(true, std::memory_order_relaxed);
                return;
            }

            for (RIImpl::State &child : children)
                pushLocal(tid, std::move(child));
        }
    }

    /**
     * Take a state off worker @a tid's own end of its queue.
     *
     * @return false if this worker has nothing left at all, which is the definition of idle that
     *         the termination protocol rests on.
     */
    bool popLocal(index tid, RIImpl::State &out) {
        Worker &worker = workers[tid];

        // The common case, and the reason for the split: my own newest state, no atomic anywhere.
        if (!worker.states.empty()) {
            out = std::move(worker.states.back());
            worker.states.pop_back();
            return true;
        }

        // Nothing private left, so take back what nobody stole. Reporting "idle" while states are
        // still sitting here would let the token finish a lap with work outstanding. The load is
        // exact rather than a hint: only the owner ever *adds* to `stealable`, so a zero it reads
        // about its own queue cannot go stale under it.
        if (worker.offered.load() == 0)
            return false;

        lockQueue(worker.busy);
        const QueueUnlock guard(worker.busy);
        if (worker.stealable.empty())
            return false;

        // From the newest end: the oldest are the biggest subtrees and are worth more to a thief.
        out = std::move(worker.stealable.back());
        worker.stealable.pop_back();
        worker.offered.store(worker.stealable.size());
        return true;
    }

    /// Push a freshly produced state onto worker @a tid's own end, publishing a batch when enough
    /// have piled up.
    void pushLocal(index tid, RIImpl::State &&state) {
        Worker &worker = workers[tid];
        worker.states.push_back(std::move(state));
        if (++worker.sinceLastPublish >= TaskGroupSize)
            coalesceIntoTask(tid);
    }

    /**
     * Make a batch of this worker's oldest states available for stealing.
     *
     * Expanding one state is far too cheap to justify the synchronization a steal costs, so states
     * become visible in batches rather than one at a time.
     */
    void coalesceIntoTask(index tid) {
        Worker &worker = workers[tid];
        worker.sinceLastPublish = 0;

        // Nothing to gain from a second group while the first is still on offer, and something to
        // lose: publishing gives away the states with the most work under them.
        if (worker.offered.load() != 0)
            return;

        // Never publish down to nothing. A worker that gives away everything goes straight back to
        // stealing, which is the cost this mechanism exists to avoid.
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
     * Try to take work from somebody else.
     *
     * Returning false is normal and simply means "nothing to steal right now"; the caller then
     * goes into the termination protocol.
     */
    bool trySteal(index thief, RIImpl::State &out) {
        // Read once and handed to pickVictim(), rather than read again in there. The draw's range
        // is `active - 2` in unsigned arithmetic, so an `active` below 2 would not give a small
        // range but a colossal one, and the index it produced would run straight off `workers`.
        // Passing the value that was just checked is what makes that impossible to reintroduce by
        // deleting a guard somewhere else.
        const count active = activeWorkers.load();
        if (active < 2)
            return false;

        for (count attempt = 0; attempt < StealAttempts; ++attempt) {
            Worker &victim = workers[pickVictim(thief, active)];

            // Two cheap rejections before the expensive one. Both are loads on a line this thief
            // never writes to.
            if (victim.offered.load() == 0)
                continue;
            if (!tryLockQueue(victim.busy))
                continue;

            const QueueUnlock guard(victim.busy);
            if (victim.stealable.empty())
                continue;

            // Somebody has work again, so any termination lap in progress is void. Set *inside*
            // the victim's critical section, which is what makes the lap argument in passToken()
            // hold: the victim can only declare itself idle after finding this queue empty, and it
            // can only find it empty by taking this lock after we have let it go.
            tokenDirty.store(true);

            out = std::move(victim.stealable.front());
            victim.stealable.pop_front();
            victim.offered.store(victim.stealable.size());
            return true;
        }

        return false;
    }

    /**
     * Choose whom to steal from: a uniformly random worker other than @a thief.
     *
     * Random choice is what keeps the load balanced without anybody having to track who is busy.
     *
     * @param thief The worker asking; never returned.
     * @param active How many workers there are. Must be at least 2 - with fewer there is no other
     *        worker to name, and the range below would underflow. @ref trySteal() is the only
     *        caller and checks exactly that before passing the value on.
     */
    index pickVictim(index thief, count active) const {
        // Thread-local already, so no seeding, no shared generator and no per-worker field.
        auto &urng = Aux::Random::getURNG();

        // Drawn from a range one short and shifted past myself, so every other worker is equally
        // likely and no draw is ever thrown away.
        std::uniform_int_distribution<index> pick(0, static_cast<index>(active) - 2);
        const index drawn = pick(urng);
        return drawn < thief ? drawn : drawn + 1;
    }

    /**
     * Move the termination token one step around the ring, on behalf of an idle worker @a tid.
     *
     * Working out that *everybody* has finished is itself a synchronization problem, and a shared
     * counter of busy workers would be contended hardest exactly when the search is winding down.
     * So a token walks a ring instead: a worker with states left simply stalls it where it is,
     * which is what makes a completed lap mean something.
     *
     * Why the lap is sound. A worker calls this only about itself and only when @ref popLocal()
     * has just reported both its queues empty, so a worker in the middle of an expansion is never
     * counted as idle. A worker with nothing can only acquire work by stealing, and every
     * successful steal sets `tokenDirty` while holding the victim's queue. Suppose a lap of
     * `activeWorkers` clean hops finished at time T and some worker j still had work then. Worker
     * j hopped during that lap, so it was empty at its own hop and must have stolen afterwards and
     * before T. That steal is ordered before the lap's final hop - it cannot be the final hop's
     * own worker, which was idle - so that hop, or an earlier one, reads the flag set and zeroes
     * the counter, and the lap never completes. Hence no lap can complete while work remains.
     *
     * `tokenDirty` and `Worker::offered` are sequentially consistent for exactly that argument:
     * they are the two things the victim looks at between the thief's store and its own hop, and
     * both are cold paths where the strong ordering costs nothing.
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

    /// Whether the search is over: a full lap of the ring with nobody having found work, or the
    /// match limit reached, or an interrupt, or an exception.
    bool quiescent() const {
        return stopped.load(std::memory_order_relaxed) || tokenHops.load() >= activeWorkers.load();
    }

    const SearchGraph *patternGraph;
    const SearchGraph *targetGraph;

    const std::vector<index> *patternNodeLabels;
    const std::vector<index> *targetNodeLabels;

    /// Computed once by the caller and read by every worker. Never modified here.
    const RIImpl::Ordering *ordering;

    /// Likewise: one shared copy rather than one per worker, which is what the RI-Ds preprocessing
    /// being lifted out of RIImpl's constructor buys. Never modified here.
    const RIImpl::Domains *domains;

    SubgraphIsomorphism::Semantics semantics;

    /// Shared by every worker. Inside the region only the non-throwing isRunning() may be
    /// called on it; ParallelRI::run() does the throwing assureRunning() once, after the join.
    Aux::SignalHandler *handler;

    /// Where the user's callback is invoked. Entered by several workers at once, and may throw;
    /// see @ref stopOnException().
    Deliver deliver;
    bool storeMatches;
    /// 0 means no limit.
    count maxMatches;
    count numWorkers;

    std::vector<Worker> workers;

    /// How many workers OpenMP actually started, which is at most `numWorkers`. Written once
    /// inside `omp single` and read by everything that has to agree on the size of the ring.
    /// Atomic rather than a plain member so that ThreadSanitizer, which cannot see through
    /// libgomp's barriers, does not report the handover as a race.
    std::atomic<count> activeWorkers;

    /// Set when the match limit is reached, the search is interrupted or a worker caught an
    /// exception, so every worker unwinds.
    std::atomic<bool> stopped;
    /// Set by the first worker that catches an exception, so exactly one worker writes `failure`.
    std::atomic<bool> failed{false};
    /// The first exception a worker caught. Written inside the region and read by run() after the
    /// join, which orders the two.
    std::exception_ptr failure;
    /// Who currently holds the termination token.
    std::atomic<index> tokenHolder;
    /// Set by any thief the moment a steal succeeds. The next hop reads and clears it, and voids
    /// the lap in progress.
    std::atomic<bool> tokenDirty;
    /// Hops since the last void. Reaching `activeWorkers` means the search is over.
    std::atomic<count> tokenHops;

    /// Running total of the matches the workers have published. Only touched when maxMatches != 0
    /// and more than one worker is running.
    std::atomic<count> published;
    /// How many matches a worker records between two publications of its local count.
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

    // Asked once, here, so that the number of workers the search actually runs and the number
    // numberOfWorkers() promised its caller cannot drift apart.
    const count numWorkers = numberOfWorkers();

    // Built once here and then shared read-only by every worker - which is the whole point of
    // hoisting the preparation out of RIImpl: one snapshot, one set of domains and one matching
    // order replace the copy each worker would otherwise build for itself. The same function serves
    // RI::run(), including the narrow refusal it carries for parallel arcs whose labels disagree,
    // which is raised before any worker exists and so leaves no thread to unwind.
    const IsomorphismDetails::RISearchSetup setup = IsomorphismDetails::prepareRISearch(
        *pattern, *target, patternNodeLabels, targetNodeLabels, patternEdgeLabels, targetEdgeLabels,
        variant, "ParallelRI");

    // The lambda is what gets the workers at the user's callback: they are not subclasses and so
    // cannot reach the protected invokeCallback() themselves, but run() can.
    ParallelRIImpl impl(
        setup.patternGraph, setup.targetGraph, patternNodeLabels, targetNodeLabels, setup.ordering,
        setup.domains, semantics, handler,
        [this](index tid, const Match &match) { return invokeCallback(tid, match); },
        storesMatches(), maxMatches, numWorkers);
    impl.run();

    // Only now, with every worker joined, is it safe to let an exception out.
    handler.assureRunning();

    finishRun(impl.takeMatches(), impl.matchesFound());
}

} // namespace NetworKit
