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
 * Abstract base class for subgraph isomorphism algorithms.
 *
 * The algorithms take a small **pattern** graph and a large **target** graph and find every
 * **match**: an injective mapping of the pattern nodes to target nodes under which every pattern
 * edge is mapped to a target edge. A match is a @ref Match indexed by pattern node, so `match[u]`
 * is the target node that pattern node @a u is mapped to.
 *
 * Two semantics are supported. Take the path `a - b - c` as the pattern and the triangle
 * `x - y - z` as the target:
 *
 * - Under @ref Semantics::MONOMORPHISM the path occurs, because every pattern edge has a
 *   corresponding target edge. The additional target edge x-z does not matter.
 * - Under @ref Semantics::INDUCED the path does not occur, because the pattern non-edge a-c must
 *   be mapped to a target non-edge as well. This is the default.
 *
 * To use an algorithm, construct it, optionally call @ref setNodeLabels(), @ref setEdgeLabels(),
 * @ref setCallback() or @ref setStoreMatches(), call @ref run(), and query the result with
 * @ref getMatches(), @ref numberOfMatches() or @ref hasMatch().
 *
 * @code
 * // How many triangles does G contain?
 * Graph triangle(3);
 * triangle.addEdge(0, 1);
 * triangle.addEdge(1, 2);
 * triangle.addEdge(2, 0);
 *
 * VF2 algo(triangle, G, SubgraphIsomorphism::Semantics::MONOMORPHISM);
 * algo.setStoreMatches(false);
 * algo.run();
 * count occurrences = algo.numberOfMatches();
 * @endcode
 *
 * This counts every triangle six times, once per automorphism of the pattern. Divide by the
 * number of automorphisms of the pattern to count unordered occurrences.
 *
 * @ref VF2 is the simple reference implementation, @ref RI is usually faster on sparse targets,
 * and @ref ParallelRI runs the RI search on several threads. They find the same set of matches,
 * but @ref ParallelRI does not guarantee any order. The search of @ref VF3 is not implemented yet.
 *
 * Notes:
 * - Edge weights are ignored.
 * - The pattern must not contain self-loops. Self-loops in the target are ignored.
 * - Parallel edges are collapsed, since the search runs on the simple graph underlying each
 *   input. See @ref setEdgeLabels() for parallel edges with different labels.
 * - The argument order is (pattern, target), which is the opposite of igraph's
 *   `igraph_subisomorphic_vf2`.
 * - Pattern and target must both be directed or both be undirected.
 * - Both graphs are held by reference and must outlive the algorithm. @ref run() rechecks the
 *   input and throws if a graph was modified in a way that invalidates it, for example if a label
 *   vector became too short.
 */
class SubgraphIsomorphism : public Algorithm {

public:
    /**
     * What counts as a match. See the class documentation for an example.
     */
    enum class Semantics : uint8_t {
        /// Pattern edges map to target edges and pattern non-edges map to target non-edges.
        INDUCED,
        /// Pattern edges map to target edges. Additional target edges are allowed.
        MONOMORPHISM,
    };

    /**
     * One match, indexed by pattern node: `match[u]` is the target node that pattern node @a u is
     * mapped to. It has `upperNodeIdBound()` entries of the pattern and holds @ref none at ids
     * that are not nodes.
     */
    using Match = std::vector<node>;

    /**
     * Callback that receives one match at a time.
     *
     * It is never called concurrently, not even by @ref ParallelRI, which serializes the calls
     * with a lock. For a parallel search the callback therefore becomes a bottleneck; use
     * @ref ParallelMatchCallback to avoid this.
     *
     * The match refers to an internal buffer that is reused for the next match. Copy it to keep
     * it.
     */
    using MatchCallback = std::function<void(const Match &)>;

    /**
     * Callback that receives one match at a time together with the id of the worker that found
     * it.
     *
     * It may be called concurrently and must therefore be thread-safe. The worker id lies in
     * `[0, numberOfWorkers())`, so per-worker accumulators need no locking:
     *
     * @code
     * std::vector<Acc> perThread(algo.numberOfWorkers());
     * algo.setCallback([&](index tid, const Match &match) { perThread[tid].add(match); });
     * @endcode
     *
     * The sequential algorithms always pass worker id 0. The match refers to an internal buffer
     * that is reused for the next match.
     */
    using ParallelMatchCallback = std::function<void(index, const Match &)>;

    ~SubgraphIsomorphism() override = default;

    /**
     * Runs the search.
     *
     * Throws `Aux::SignalHandler::InterruptException` if the search is interrupted with CTRL+C. If
     * a callback throws, the search stops and the exception is rethrown; a parallel search first
     * stops all workers. In both cases the algorithm is not finished afterwards, so
     * @ref getMatches(), @ref numberOfMatches() and @ref hasMatch() throw. Matches that were
     * already passed to a callback remain delivered.
     */
    void run() override = 0;

    /**
     * Returns the number of workers @ref run() uses. The worker id passed to a
     * @ref ParallelMatchCallback is smaller than this number. A parallel algorithm reads the
     * global thread count, so call this after `Aux::setNumberOfThreads()`.
     *
     * @return the number of workers; 1 for the sequential algorithms.
     */
    virtual count numberOfWorkers() const { return 1; }

    /**
     * Restricts matches to map every pattern node to a target node with the same label. The label
     * @ref none is a wildcard that matches any label. Passing two empty vectors removes the node
     * labels. Node labels and edge labels are independent of each other.
     *
     * Call this before @ref run().
     *
     * @param patternNodeLabels Labels of the pattern nodes, indexed by node id. Needs at least
     * `upperNodeIdBound()` entries of the pattern.
     * @param targetNodeLabels Labels of the target nodes, indexed by node id. Needs at least
     * `upperNodeIdBound()` entries of the target.
     * @throws std::runtime_error if a vector is too short.
     */
    void setNodeLabels(const std::vector<index> &patternNodeLabels,
                       const std::vector<index> &targetNodeLabels);

    /**
     * Restricts matches to map every pattern edge to a target edge with the same label. The label
     * @ref none is a wildcard that matches any label. Passing two empty vectors removes the edge
     * labels. Node labels and edge labels are independent of each other.
     *
     * Both graphs need edge ids, see `Graph::indexEdges()`. Parallel edges with different labels
     * cannot be collapsed into one edge, so @ref run() throws `std::runtime_error` for such input.
     * @ref VF3 does not support edge labels and throws from @ref run() if they are set.
     *
     * Call this before @ref run().
     *
     * @param patternEdgeLabels Labels of the pattern edges, indexed by edge id. Needs at least
     * `upperEdgeIdBound()` entries of the pattern.
     * @param targetEdgeLabels Labels of the target edges, indexed by edge id. Needs at least
     * `upperEdgeIdBound()` entries of the target.
     * @throws std::runtime_error if a graph has no edge ids or a vector is too short.
     */
    void setEdgeLabels(const std::vector<index> &patternEdgeLabels,
                       const std::vector<index> &targetEdgeLabels);

    /**
     * Passes every match to @a callback instead of storing it. @ref getMatches() then throws,
     * while @ref numberOfMatches() and @ref hasMatch() keep working. Replaces a callback that was
     * set earlier, of either form.
     *
     * Call this before @ref run().
     *
     * @param callback Called once per match, never concurrently. See @ref MatchCallback.
     */
    void setCallback(MatchCallback callback);

    /**
     * Like @ref setCallback(MatchCallback), but for a thread-safe callback that also receives the
     * worker id and may be called concurrently. This lets @ref ParallelRI use all its workers.
     *
     * Call this before @ref run().
     *
     * @param callback Called once per match; must be thread-safe. See @ref ParallelMatchCallback.
     */
    void setCallback(ParallelMatchCallback callback);

    /**
     * Sets whether matches are stored. Pass false to only count them. @ref getMatches() then
     * throws, while @ref numberOfMatches() and @ref hasMatch() keep working. Has no effect if a
     * callback is set, since matches are then never stored.
     *
     * Call this before @ref run().
     *
     * @param storeMatches Whether to store matches for @ref getMatches(). Default: true.
     */
    void setStoreMatches(bool storeMatches);

    /**
     * Returns all matches found by @ref run().
     *
     * Throws `std::runtime_error` if @ref run() has not finished, if a callback was set, or if
     * @ref setStoreMatches(false) was called.
     *
     * @return the matches, each indexed by pattern node.
     */
    const std::vector<Match> &getMatches() const;

    /**
     * Returns the number of matches found, whether or not they were stored.
     *
     * With a cap on the number of matches, the result is at most the cap, with one exception: a
     * parallel search with a callback may deliver a few matches beyond the cap before all workers
     * stop, and these matches are counted too.
     *
     * @return the number of matches found.
     */
    count numberOfMatches() const;

    /**
     * Returns whether at least one match was found. To answer only this question, construct the
     * algorithm with `maxMatches = 1`, so that the search stops at the first match.
     *
     * @return true if at least one match was found.
     */
    bool hasMatch() const;

protected:
    // A sequential run() calls prepareRun(), then reportMatch() for every match until it returns
    // false, and finally finishRun(). A parallel run() calls prepareRun() and passes every match to
    // invokeCallback(). If invokeCallback() returns false and storesMatches() is true, the worker
    // stores the match in its own buffer. After the workers have joined, run() calls
    // finishRun(matches, found). Inside a parallel region, only Aux::SignalHandler::isRunning() may
    // be polled; assureRunning() is called after the join. ParallelRI.cpp is the worked example.

    /**
     * @param pattern The pattern graph to look for.
     * @param target The target graph to look in.
     * @param semantics Whether matches must be induced.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     */
    SubgraphIsomorphism(const Graph &pattern, const Graph &target, Semantics semantics,
                        count maxMatches);

    /**
     * @return true if node labels were set with @ref setNodeLabels().
     */
    bool isNodeLabelled() const noexcept { return !patternNodeLabels.empty(); }

    /**
     * @return true if edge labels were set with @ref setEdgeLabels(). An algorithm that does not
     * support edge labels must throw from @ref run() in this case.
     */
    bool isEdgeLabelled() const noexcept { return !patternEdgeLabels.empty(); }

    /**
     * @return true if a callback of either form was set.
     */
    bool hasCallback() const noexcept {
        return static_cast<bool>(callback) || static_cast<bool>(parallelCallback);
    }

    /**
     * @return true if a @ref MatchCallback was set, which must not be called concurrently.
     */
    bool hasSerialCallback() const noexcept { return static_cast<bool>(callback); }

    /**
     * @return true if found matches have to be stored, that is, no callback is set and
     * @ref setStoreMatches(false) was not called.
     */
    bool storesMatches() const noexcept { return storeMatches && !hasCallback(); }

    /**
     * Resets the results of an earlier run and rechecks the input, since the graphs may have
     * changed since construction. Throws `std::runtime_error` if the input is no longer valid.
     * Call this at the start of @ref run().
     */
    void prepareRun();

    /**
     * Records one match of a sequential search: passes it to the callback if one was set, stores
     * it otherwise, and counts it. Not thread-safe; a parallel search uses @ref invokeCallback()
     * instead.
     *
     * @param match The mapping, indexed by pattern node.
     * @return false once the cap on the number of matches is reached, in which case the search
     * must stop; true otherwise.
     */
    bool reportMatch(const Match &match);

    /**
     * @return how many more matches a sequential search may report, or @ref none if there is no
     * cap. Use it to size a buffer; use the return value of @ref reportMatch() to decide when to
     * stop.
     */
    count remainingMatches() const noexcept {
        if (maxMatches == 0)
            return none;
        return matchCount >= maxMatches ? 0 : maxMatches - matchCount;
    }

    /**
     * Passes one match to the callback. May be called from any thread: a
     * @ref ParallelMatchCallback is called directly, a @ref MatchCallback under a lock. Does not
     * count or store the match.
     *
     * @param tid Worker id in `[0, numberOfWorkers())`.
     * @param match The mapping, indexed by pattern node.
     * @return true if a callback received the match, in which case the caller must not store it.
     */
    bool invokeCallback(index tid, const Match &match);

    /**
     * Marks a sequential run as finished. Call this at the end of @ref run().
     */
    void finishRun();

    /**
     * Marks a parallel run as finished. Adopts the merged matches of all workers and, unless a
     * callback is set, trims them to the cap on the number of matches. Call this once, after all
     * workers have joined.
     *
     * @param matches The concatenated matches of all workers. Empty if no matches were stored.
     * @param found The total number of matches found, stored or not.
     */
    void finishRun(std::vector<Match> &&matches, count found);

    /// The graph to look for.
    const Graph *pattern;
    /// The graph to look in.
    const Graph *target;

    /// Labels of the pattern nodes, indexed by node id. Empty unless setNodeLabels() was called.
    std::vector<index> patternNodeLabels;
    /// Labels of the target nodes, indexed by node id. Empty unless setNodeLabels() was called.
    std::vector<index> targetNodeLabels;

    /// Labels of the pattern edges, indexed by edge id. Empty unless setEdgeLabels() was called.
    std::vector<index> patternEdgeLabels;
    /// Labels of the target edges, indexed by edge id. Empty unless setEdgeLabels() was called.
    std::vector<index> targetEdgeLabels;

    Semantics semantics;
    /// Maximum number of matches; 0 means no limit.
    count maxMatches;

private:
    /**
     * Throws `std::runtime_error` if pattern and target differ in directedness or the pattern has
     * self-loops. Called from the constructor and from @ref prepareRun().
     */
    void validateInput() const;

    /**
     * Throws `std::runtime_error` if the node label vectors are too short for the graphs. Two
     * empty vectors are always accepted.
     */
    void validateNodeLabels(const std::vector<index> &patternNodeLabels,
                            const std::vector<index> &targetNodeLabels) const;

    /**
     * Throws `std::runtime_error` if a graph has no edge ids or the edge label vectors are too
     * short for the graphs. Two empty vectors are always accepted.
     */
    void validateEdgeLabels(const std::vector<index> &patternEdgeLabels,
                            const std::vector<index> &targetEdgeLabels) const;

    std::vector<Match> result;

    MatchCallback callback;
    ParallelMatchCallback parallelCallback;

    /// Serializes the calls of `callback` in invokeCallback().
    std::mutex reportMutex;

    count matchCount;
    bool storeMatches;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_
