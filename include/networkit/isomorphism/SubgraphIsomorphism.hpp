#ifndef NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_
#define NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_

#include <cstdint>
#include <functional>
#include <mutex>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup isomorphism
 * Common base class of the subgraph isomorphism algorithms. Start reading here.
 *
 * ## What these algorithms do
 *
 * You have a small graph, called the **pattern**, and a big graph, called the **target**. You
 * want to find every place in the target where the pattern occurs. Each such place is called a
 * **match**.
 *
 * The classic example is counting triangles. The pattern is a triangle: three nodes 0, 1, 2 with
 * the edges 0-1, 1-2, 2-0. The target is your network. Every match tells you three target nodes
 * that form a triangle.
 *
 * A match is a @ref Match - a `std::vector<node>` **indexed by pattern node**: `match[u]` is the
 * target node that pattern node @a u was mapped to. So for the triangle pattern above, a match of
 * `{17, 42, 8}` means pattern node 0 sits on target node 17, pattern node 1 on target node 42 and
 * pattern node 2 on target node 8 - and the target really does contain the edges 17-42, 42-8 and
 * 8-17.
 *
 * Two pattern nodes are never mapped to the same target node. In mathematical terms the mapping
 * is injective.
 *
 * ## Induced or not: the one decision you must get right
 *
 * Consider a pattern that is a path of three nodes, `a - b - c`, and a target that is a triangle
 * `x - y - z` (so all three target edges exist).
 *
 * - Under @ref Semantics::MONOMORPHISM the answer is "yes, the path occurs": every pattern edge
 *   (a-b and b-c) has a corresponding target edge. That the target additionally has the edge x-z
 *   does not matter.
 * - Under @ref Semantics::INDUCED the answer is "no": the pattern has **no** edge between a and
 *   c, so the target must have **no** edge between the nodes they map to. But x-z exists, so this
 *   is not an induced occurrence.
 *
 * Rule of thumb: if your question is "does this shape occur somewhere", you want MONOMORPHISM.
 * If your question is "does exactly this shape occur, with nothing extra", you want INDUCED.
 * Path and tree queries are almost always MONOMORPHISM questions. The default is INDUCED because
 * that is the stricter, less surprising answer.
 *
 * ## How to use it
 *
 * Four steps, always in this order:
 *
 * 1. **Construct** with the pattern, the target, the semantics and an optional cap on how many
 *    matches you want.
 * 2. **Configure**, if you need to: @ref setNodeLabels() to restrict matches to nodes that carry
 * the same label, @ref setEdgeLabels() to do the same for edges, @ref setCallback() to be handed
 *    each match as it is found instead of collecting them, @ref setStoreMatches() to only count
 *    matches without keeping them.
 * 3. **Run** with @ref run().
 * 4. **Query** with @ref getMatches(), @ref numberOfMatches() or @ref hasMatch().
 *
 * @code
 * // How many triangles does G contain?
 * Graph triangle(3);
 * triangle.addEdge(0, 1);
 * triangle.addEdge(1, 2);
 * triangle.addEdge(2, 0);
 *
 * VF2 algo(triangle, G, SubgraphIsomorphism::Semantics::MONOMORPHISM);
 * algo.setStoreMatches(false);   // we only want the count
 * algo.run();
 * count occurrences = algo.numberOfMatches();
 * @endcode
 *
 * Note that this counts every triangle six times, once per way of labelling its corners. That is
 * inherent to subgraph matching, not a quirk of this implementation: divide by the number of
 * automorphisms of your pattern if you want unordered occurrences.
 *
 * ## Which algorithm should I use?
 *
 * | Class          | Use it when                                                              |
 * | -------------- | ------------------------------------------------------------------------ |
 * | @ref VF2       | Small inputs, or as the reference to check the others against.            |
 * | @ref VF3       | **Not implemented yet** - @ref VF3::run() throws. See its documentation.  |
 * | @ref RI        | Sparse targets and biological networks; usually the fastest sequentially. |
 * | @ref ParallelRI| Same as RI, but you have cores to spare and the search is long.           |
 *
 * Every algorithm that is implemented returns exactly the same set of matches. They differ only in
 * how fast they get there, and @ref ParallelRI additionally does not promise any particular order.
 * @ref VF3 is the exception, and only because its search is still unwritten: the class, the
 * documentation and the input validation are in place, but @ref VF3::run() throws a
 * `std::logic_error` rather than searching. Nothing here is planned differently for it - when it
 * lands it will answer exactly what the other three do.
 *
 * ## Things that quietly do not apply
 *
 * - **Edge weights are ignored.** Only the structure of the two graphs is matched.
 * - **Self-loops in the pattern are rejected** by the constructor. Self-loops in the *target* are
 *   fine and are simply never used, because both semantics only constrain pairs of distinct
 *   nodes.
 * - **Parallel edges are collapsed.** The search runs on the simple graph underlying each input.
 *   Multiplicity cannot change which matches exist - both semantics only ever ask whether a pair
 *   of distinct nodes is joined, which is a yes-or-no question that a second copy of an edge does
 *   not change - so this normalization is what keeps degrees honest and stops the same candidate
 *   being enumerated, and reported, twice. `Graph::addEdge()` permits parallel edges by default,
 *   so this is reachable without doing anything unusual. Edge labels are the one exception: two
 *   parallel edges carrying *different* labels cannot collapse into one arc that stands for both,
 *   and @ref setEdgeLabels() documents how that is refused rather than papered over.
 * - **The argument order is (pattern, target)**, which is the opposite of igraph's
 *   `igraph_subisomorphic_vf2`. Swapping the two compiles fine and answers a different question.
 * - Pattern and target must agree on directedness; mixing them throws from the constructor.
 * - **The two graphs are held by reference, not copied**, so both must outlive the algorithm
 *   object. They may be modified in between, but everything configured beforehand is rechecked
 *   when @ref run() starts: adding a node or an edge after @ref setNodeLabels() or
 *   @ref setEdgeLabels() leaves the label vector too short, and @ref run() then throws instead of
 *   reading past its end. Set the labels again at the new sizes and the run succeeds.
 *
 * ## What this class contributes
 *
 * Everything that does not depend on *how* the search is done: the two graphs, the optional
 * labels, the semantics, the result store, the callbacks, and the cap on the number of matches.
 * A concrete algorithm therefore only has to implement @ref run(). See the protected section for
 * the protocol that @ref run() has to follow.
 *
 * ## Ideas for later
 *
 * Things left out on purpose, recorded so the next reader can tell a decision from an
 * oversight.
 *
 * - **Parallel edges whose labels disagree.** @ref setEdgeLabels() refuses these today. The
 *   search runs on a snapshot that collapses repeated neighbours, which is what keeps degrees
 *   honest and stops the same candidate being enumerated twice; collapsing two arcs whose
 *   labels differ would quietly answer a different question. Supporting them does *not* mean
 *   giving the collapse up. It means giving each collapsed arc a **set** of labels instead of
 *   one, and reading "pattern edge with label l matches" as "the image pair carries some
 *   compatible label". The structural snapshot is untouched, so degree pruning and candidate
 *   enumeration are unaffected and no match can be reported twice; only the label comparison
 *   widens from equality to set intersection. That is enough for multi-relational graphs, where
 *   a pair of nodes is joined by several differently-typed edges.
 *
 *   The stronger reading - every pattern edge matched by its *own distinct* target edge - is a
 *   different problem, not a bigger version of this one. It turns each accepted node pair into
 *   a small bipartite matching, and a match is an array of node images, so it could report
 *   *that* an assignment exists but never *which*. `Graph::edgeId(u, v)` cannot name the second
 *   parallel edge either. That variant needs its own result type and its own reference matcher.
 *
 * - **Labels as graph attributes.** @ref setNodeLabels() and @ref setEdgeLabels() both take flat
 *   vectors, which predates `Attributes.hpp`. Naming an `index`-typed node or edge attribute
 *   instead would be a friendlier API and would let labels travel with the graph through copies.
 *   It would change only how labels are *supplied*: the search probes them in its innermost
 *   loop, so a setter would still materialise them into flat vectors once rather than paying a
 *   hash lookup per probe. One gap to close first - `AdjListGraph`'s typed convenience wrappers
 *   cover `int`, `double` and `std::string` only, so an `index`-typed attribute is not reachable
 *   from Python until those gain an `index`-typed sibling.
 */
class SubgraphIsomorphism : public Algorithm {

public:
    /**
     * What counts as a match. See the class documentation for a worked example of the
     * difference; the short version is below.
     */
    enum class Semantics : uint8_t {
        /// Pattern edges must map to target edges **and** pattern non-edges to target non-edges.
        INDUCED,
        /// Only pattern edges must map to target edges. Extra target edges are allowed.
        MONOMORPHISM,
    };

    /**
     * One match, in the shape this module defines everywhere.
     *
     * Indexed by **pattern node**: `match[u]` is the target node that pattern node @a u was
     * mapped to. Sized to the pattern's `upperNodeIdBound()`, holding @ref none at ids that are
     * not nodes, so a pattern with removed ids still indexes directly.
     *
     * A plain alias rather than a distinct type, so any `std::vector<node>` still passes.
     */
    using Match = std::vector<node>;

    /**
     * Callback handed one match at a time, instead of collecting them all in memory.
     *
     * The argument is the mapping indexed by pattern node, exactly what @ref getMatches() would
     * have stored.
     *
     * This form is **never called from two threads at once**, not even by @ref ParallelRI, which
     * serializes it behind a lock; it needs no locking of its own. The price is that it becomes a
     * bottleneck: a parallel search cannot go faster than this callback runs. If that matters, use
     * @ref ParallelMatchCallback instead.
     *
     * The reference points at an internal buffer that is reused for the next match. Copy the
     * vector if you need to keep it.
     */
    using MatchCallback = std::function<void(const Match &)>;

    /**
     * Same as @ref MatchCallback, but it also receives the id of the worker that found the match.
     *
     * This form **may be called from several threads at once** and must therefore be thread-safe.
     * The worker id is in `[0, numberOfWorkers)`, which lets you give every worker its own slot
     * and avoid locking entirely:
     *
     * @code
     * std::vector<Acc> perThread(algo.numberOfWorkers());
     * algo.setCallback([&](index tid, const Match &match) {
     *     perThread[tid].add(match);   // no lock needed, each tid owns its slot
     * });
     * @endcode
     *
     * The sequential algorithms always pass `tid == 0`, so code written against this form works
     * unchanged with @ref VF2, @ref VF3 and @ref RI.
     *
     * The reference points at an internal buffer that is reused for the next match.
     */
    using ParallelMatchCallback = std::function<void(index, const Match &)>;

    ~SubgraphIsomorphism() override = default;

    /**
     * Run the search. Implemented by each concrete algorithm.
     *
     * Throws `Aux::SignalHandler::InterruptException` if the search was interrupted with CTRL+C.
     * The algorithm is then left **not finished**: @ref Algorithm::hasFinished() stays false and
     * @ref getMatches(), @ref numberOfMatches() and @ref hasMatch() all throw. Matches already
     * handed to a callback stay handed over - a search cannot take them back. This is what every
     * interruptible algorithm in NetworKit does, so nothing here is a special case.
     *
     * An exception thrown by a callback ends the run the same way: the search stops, the algorithm
     * is left not finished, and the exception comes out of run(). A parallel search first stops
     * every worker and then rethrows the first exception any of them caught.
     */
    void run() override = 0;

    /**
     * How many workers @ref run() will use.
     *
     * This is the bound the worker id handed to a @ref ParallelMatchCallback stays below, so it is
     * what you size a per-worker accumulator by. The sequential algorithms return 1.
     *
     * Ask *after* any `Aux::setNumberOfThreads()` call: a parallel algorithm reads the thread
     * count when @ref run() starts, so changing the setting in between changes the answer.
     */
    virtual count numberOfWorkers() const { return 1; }

    /**
     * Only accept matches that map like-labelled nodes onto each other.
     *
     * This is the node-side of the two label APIs; @ref setEdgeLabels() is its edge-side
     * counterpart, and the two are independent - setting one does not disturb the other.
     *
     * Both vectors are indexed by node id, so `patternNodeLabels[u]` is the label of pattern node
     * @a u. A pattern node may only be mapped to a target node carrying the same label. The
     * special value @ref none acts as a wildcard and matches any label.
     *
     * Both vectors must have at least `upperNodeIdBound()` entries for their graph, otherwise
     * `std::runtime_error` is thrown. Passing two empty vectors clears the node labels and puts
     * the algorithm back into the node-unlabelled search.
     *
     * Node labels usually make the search much faster, because most candidate pairs can be
     * rejected without looking at the graph structure at all. @ref VF3 gets the most out of them.
     *
     * If you already have a `Partition` - from community detection, say - pass its
     * `getVector()`.
     *
     * Call this before @ref run().
     *
     * @param patternNodeLabels Labels of the pattern nodes, indexed by node id.
     * @param targetNodeLabels Labels of the target nodes, indexed by node id.
     */
    void setNodeLabels(const std::vector<index> &patternNodeLabels,
                       const std::vector<index> &targetNodeLabels);

    /**
     * Only accept matches that map like-labelled edges onto each other.
     *
     * This is the edge-side counterpart of @ref setNodeLabels(), and it is what a multi-relational
     * graph needs: when the same pair of nodes can be joined by a "cites", a "co-authors" and a
     * "rebuts" edge, "this pattern edge must land on a target edge of the same kind" is not
     * expressible with node labels alone.
     *
     * Both vectors are indexed by **edge id**, not by node id, so both graphs must have had
     * `indexEdges()` called on them and each vector must have at least `upperEdgeIdBound()`
     * entries for its graph; otherwise `std::runtime_error` is thrown. The special value
     * @ref none acts as a wildcard and matches any label, on either side. Passing two empty
     * vectors clears the edge labels, exactly as @ref setNodeLabels() does for node labels.
     *
     * Two limits, both deliberate and both raised as `std::runtime_error` rather than answered
     * wrongly. **Parallel edges whose labels disagree are refused** - the search runs on a
     * snapshot that collapses repeated neighbours, so there is no one label for the collapsed
     * arc to carry; this is detected while that snapshot is built, so it surfaces from
     * @ref run() rather than from here. And **not every algorithm understands edge labels yet**:
     * @ref VF2 and @ref VF3 refuse them outright from @ref run(). See "Ideas for later" in the
     * class documentation for what lifting the first limit would take.
     *
     * Call this before @ref run().
     *
     * @param patternEdgeLabels Labels of the pattern edges, indexed by edge id.
     * @param targetEdgeLabels Labels of the target edges, indexed by edge id.
     */
    void setEdgeLabels(const std::vector<index> &patternEdgeLabels,
                       const std::vector<index> &targetEdgeLabels);

    /**
     * Hand each match to @a callback as it is found, rather than collecting them.
     *
     * Use this when there may be too many matches to hold in memory. @ref getMatches() throws
     * afterwards, but @ref numberOfMatches() and @ref hasMatch() keep working. Replaces any
     * callback set earlier, of either form.
     *
     * This callback is never invoked concurrently. See @ref MatchCallback for what that costs.
     *
     * Call this before @ref run().
     *
     * @param callback Called once per match.
     */
    void setCallback(MatchCallback callback);

    /**
     * Like the other @ref setCallback(), but for a thread-safe callback that also receives the
     * worker id and may be called from several threads at once.
     *
     * This is the form that lets @ref ParallelRI use all its workers. Replaces any callback set
     * earlier, of either form.
     *
     * Call this before @ref run().
     *
     * @param callback Called once per match; must be thread-safe.
     */
    void setCallback(ParallelMatchCallback callback);

    /**
     * Choose whether the matches are kept at all.
     *
     * Pass false when you only want to know *how many* matches there are, as in motif or graphlet
     * counting. Nothing is allocated per match and there is no merge pass at the end, which makes
     * it the cheapest way to ask the question. @ref getMatches() then throws, while
     * @ref numberOfMatches() and @ref hasMatch() keep working.
     *
     * Has no effect if a callback is set, since a callback already means nothing is stored.
     *
     * Call this before @ref run().
     *
     * @param storeMatches Whether to keep matches for @ref getMatches(). Default: true.
     */
    void setStoreMatches(bool storeMatches);

    /**
     * All matches found, each one a vector indexed by pattern node.
     *
     * Throws `std::runtime_error` if the matches were never stored, which happens when a callback
     * was set or @ref setStoreMatches(false) was called, and if @ref run() has not been called
     * yet.
     *
     * @return the matches, each a vector of target node ids.
     */
    const std::vector<Match> &getMatches() const;

    /**
     * How many matches were found. Works regardless of whether they were stored.
     *
     * If a cap was requested, this is capped too and is therefore not the total number of matches
     * in the graph. The cap is exact, with one exception: a parallel search that streams to a
     * callback may have handed over a small bounded number of extra matches before every worker
     * noticed that the cap had been reached, and this then **counts those too**, so it can come
     * back slightly above the cap. Nothing is trimmed on that path on purpose - the matches really
     * were delivered, and a count that hid them would disagree with what the callback saw. Every
     * other combination is exact.
     *
     * @return the number of matches reported.
     */
    count numberOfMatches() const;

    /**
     * Whether at least one match was found.
     *
     * If all you want is a yes/no answer, construct the algorithm with `maxMatches = 1` so the
     * search stops at the first match instead of enumerating all of them.
     *
     * @return true if there was at least one match.
     */
    bool hasMatch() const;

protected:
    // ---------------------------------------------------------------------------------------
    // Everything below is for people implementing a new algorithm in this module.
    //
    // The protocol a *sequential* run() must follow - VF2, VF3 and RI all look like this:
    //
    //   1. Call prepareRun(), which throws away the results of any earlier run.
    //   2. Search. Whenever a complete mapping is found, call reportMatch(). Stop as soon as it
    //      returns false; that is how the cap on the number of matches is honoured.
    //   3. Call finishRun(). It applies the cap and marks the run finished.
    //
    // A *parallel* run() cannot use reportMatch(), which counts and stores without any locking
    // and is therefore single-threaded only. It instead does its own accumulation, which it
    // needs per-worker state for anyway, and hands the merged output over at the end:
    //
    //   1. Call prepareRun().
    //   2. Search. On a complete mapping call invokeCallback(tid, match) - that part *is*
    //      thread-safe - and, if it returns false and storesMatches() says so, append the match
    //      to that worker's own buffer. Count and coordinate the cap privately.
    //   3. After the workers have joined, call finishRun(mergedMatches, totalFound).
    //
    // ParallelRI.cpp is the worked example: its workers already own a padded slot each, so the
    // buffers and counters go there rather than into this class, where three sequential
    // algorithms would pay for machinery they never use.
    //
    // Interruption works exactly as it does everywhere else in NetworKit, with no machinery of
    // this module's own: hold an `Aux::SignalHandler` local to run() and call assureRunning() on
    // it, as MaximalCliques does. A *parallel* search must instead call the non-throwing
    // isRunning() inside the region and assureRunning() once after the workers join, because an
    // exception escaping an OpenMP structured block is undefined behaviour - see Betweenness.cpp,
    // which is the pattern to copy.
    //
    // The same rule covers an exception thrown by a callback, which a parallel search reaches
    // from inside the region. Catch it there, stop the other workers, and rethrow it after the
    // join. ParallelRIImpl::stopOnException() in ParallelRI.cpp does exactly that.
    //
    // Use isNodeLabelled() to find out whether node labels are in play, and isEdgeLabelled() for
    // edge labels. The module's rule about the latter is: refuse edge labels outright in the
    // algorithms that will not understand them, and refuse only disagreeing parallel edges in the
    // ones that will. VF2 and VF3 do the former, RI and ParallelRI the latter.
    //
    // The search itself belongs in a separate implementation class in the .cpp file, so that the
    // public header never grows a member. Those classes are not subclasses and so cannot reach
    // anything here: give them an IsomorphismDetails::MatchReporter built from a lambda in
    // run(), which can. See VF2.cpp for the pattern.
    // ---------------------------------------------------------------------------------------

    /**
     * @param pattern The pattern graph to look for.
     * @param target The target graph to look in.
     * @param semantics Whether matches must be induced.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     */
    SubgraphIsomorphism(const Graph &pattern, const Graph &target, Semantics semantics,
                        count maxMatches);

    /**
     * @return true if @ref setNodeLabels() was used, so the search has to compare node labels.
     *
     * Independent of @ref isEdgeLabelled(): either, both or neither may be true.
     */
    bool isNodeLabelled() const noexcept { return !patternNodeLabels.empty(); }

    /**
     * @return true if @ref setEdgeLabels() was used, so the search has to compare edge labels.
     *
     * An algorithm that cannot honour edge labels must throw from @ref run() when this is true.
     * Reporting matches that violate an edge label the caller asked for, with nothing to say so,
     * is the one failure worse than refusing.
     */
    bool isEdgeLabelled() const noexcept { return !patternEdgeLabels.empty(); }

    /**
     * @return true if a callback of either form was set.
     */
    bool hasCallback() const noexcept {
        return static_cast<bool>(callback) || static_cast<bool>(parallelCallback);
    }

    /**
     * @return true if the user's callback must not be called from two threads at once. Handled
     * by @ref invokeCallback(); exposed for algorithms that want to decide something else on it.
     */
    bool hasSerialCallback() const noexcept { return static_cast<bool>(callback); }

    /**
     * @return true if a search has to keep the matches it finds, rather than only count them.
     * False when a callback is set, since a callback already means nothing is stored.
     */
    bool storesMatches() const noexcept { return storeMatches && !hasCallback(); }

    /**
     * Step 1 of the protocol: forget the results of any earlier run, then recheck the input.
     *
     * The recheck matters because both graphs are held by pointer: anything the caller changed
     * since construction is caught here, throwing `std::runtime_error`, rather than being read
     * out of bounds mid-search. It runs after the reset, so a rejected run leaves nothing
     * queryable.
     */
    void prepareRun();

    /**
     * Record one match, for a sequential search.
     *
     * Hands the match to the callback if one was set, stores it otherwise, and counts it.
     *
     * Does no locking and touches plain members, so it must not be called from more than one
     * thread. A parallel search uses @ref invokeCallback() and its own buffers instead.
     *
     * @param match The mapping, indexed by pattern node. Copied if it is being stored.
     * @return true to keep searching, false once the requested number of matches is reached. A
     *         search must stop as soon as this returns false.
     */
    bool reportMatch(const Match &match);

    /**
     * How many more matches a sequential search may report, or @ref none if there is no cap.
     *
     * Use it to size a `reserve()`; do **not** use it to decide when to stop, that is what the
     * return value of @ref reportMatch() is for.
     */
    count remainingMatches() const noexcept {
        if (maxMatches == 0)
            return none;
        return matchCount >= maxMatches ? 0 : maxMatches - matchCount;
    }

    /**
     * Hand one match to the user's callback, from any thread.
     *
     * This is the part of reporting a parallel search cannot do for itself, because a
     * @ref MatchCallback promises never to be entered twice at once and the reporting call sits
     * deep inside the search, far from the code that knows how many threads are running. A
     * @ref ParallelMatchCallback is invoked directly with @a tid; a @ref MatchCallback is invoked
     * under a lock.
     *
     * Counting and storing are deliberately *not* done here: a parallel search keeps those
     * per-worker, so that the common case touches no memory another thread writes to.
     *
     * @param tid Worker id in `[0, numberOfWorkers())`.
     * @param match The mapping, indexed by pattern node.
     * @return true if a callback consumed the match, in which case the caller must not store it.
     */
    bool invokeCallback(index tid, const Match &match);

    /**
     * Step 3 of the protocol, sequential: apply the cap and mark the run as finished.
     */
    void finishRun();

    /**
     * Step 3 of the protocol, parallel: adopt the merged per-worker output, apply the cap, mark
     * the run as finished.
     *
     * Call once, after every worker has joined.
     *
     * @param matches The workers' buffers concatenated. Empty when nothing was stored.
     * @param found How many matches were reported in total, stored or not.
     */
    void finishRun(std::vector<Match> &&matches, count found);

    /// The graph we are looking for. Never null, never modified.
    const Graph *pattern;
    /// The graph we are looking in. Never null, never modified.
    const Graph *target;

    /// Empty unless setNodeLabels() was used. Indexed by node id.
    std::vector<index> patternNodeLabels;
    /// Empty unless setNodeLabels() was used. Indexed by node id.
    std::vector<index> targetNodeLabels;

    /// Empty unless setEdgeLabels() was used. Indexed by **edge** id, not node id.
    std::vector<index> patternEdgeLabels;
    /// Empty unless setEdgeLabels() was used. Indexed by **edge** id, not node id.
    std::vector<index> targetEdgeLabels;

    Semantics semantics;
    /// 0 means no limit.
    count maxMatches;

private:
    /**
     * Reject inputs no search can handle, throwing `std::runtime_error`.
     *
     * Called from the constructor, so bad input fails immediately rather than deep inside
     * @ref run(), and again from @ref prepareRun(), because the graphs are held by pointer and
     * the caller may have changed them in between. Rejects a directedness mismatch and
     * self-loops in the pattern.
     *
     * Deliberately not virtual: it first runs while the derived part of the object is still
     * uninitialized, so an override would never be reached.
     */
    void validateInput() const;

    /**
     * Reject node label vectors that do not fit the graphs, throwing `std::runtime_error`.
     *
     * Takes the vectors as arguments rather than reading the members, so @ref setNodeLabels()
     * can check *before* assigning and a rejected call leaves the object untouched, while
     * @ref prepareRun() passes the members to recheck them against the current graphs. Two empty
     * vectors are the documented way to clear the labels and are always accepted.
     */
    void validateNodeLabels(const std::vector<index> &patternNodeLabels,
                            const std::vector<index> &targetNodeLabels) const;

    /**
     * Reject edge label vectors that do not fit the graphs, throwing `std::runtime_error`.
     *
     * The edge-side counterpart of @ref validateNodeLabels(), and used the same way. Also
     * rejects graphs without edge ids, since the vectors are indexed by edge id.
     *
     * Partly belt and braces: the algorithms that honour edge labels hand the vector to
     * `SearchGraph`, whose constructor checks it again while building the snapshot. The node
     * label vector has no such second reader, which is why the recheck matters more there.
     */
    void validateEdgeLabels(const std::vector<index> &patternEdgeLabels,
                            const std::vector<index> &targetEdgeLabels) const;

    std::vector<Match> result;

    MatchCallback callback;
    ParallelMatchCallback parallelCallback;

    /// Guards `callback` in invokeCallback(), which only a parallel search calls.
    std::mutex reportMutex;

    count matchCount;
    bool storeMatches;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_
