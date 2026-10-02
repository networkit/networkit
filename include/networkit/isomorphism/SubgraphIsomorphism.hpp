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
 * The algorithms find every match of a small pattern graph in a large target graph. A match maps
 * the pattern nodes injectively to target nodes and every pattern edge to a target edge. Under
 * @ref Semantics::INDUCED, the default, it also maps every pattern non-edge to a target non-edge.
 * The search finds an occurrence once per automorphism of the pattern, so a triangle occurs six
 * times.
 *
 * @ref VF2 is the reference implementation, and @ref RI is usually faster on sparse targets.
 * @ref ParallelRI runs RI on several threads and reports the matches in no fixed order.
 *
 * Pattern and target must both be directed or both be undirected, and the pattern must not contain
 * self-loops. The search ignores edge weights and target self-loops, and it collapses parallel
 * edges. The algorithm holds both graphs by reference, so they must outlive it.
 */
class SubgraphIsomorphism : public Algorithm {

public:
    enum class Semantics : uint8_t {
        /// Pattern edges map to target edges and pattern non-edges map to target non-edges.
        INDUCED,
        /// Pattern edges map to target edges. Additional target edges are allowed.
        MONOMORPHISM,
    };

    /// One match: `match[u]` is the image of pattern node u, or @ref none if u is not a node.
    using Match = std::vector<node>;

    /**
     * Receives one match at a time, never concurrently. The match refers to an internal buffer
     * that the next match overwrites. @ref ParallelRI serializes the calls with a lock, which a
     * @ref ParallelMatchCallback avoids.
     */
    using MatchCallback = std::function<void(const Match &)>;

    /**
     * Receives one match at a time together with the id of the worker that found it, in
     * `[0, numberOfWorkers())`. It may be called concurrently and must be thread-safe. The
     * sequential algorithms pass worker id 0.
     */
    using ParallelMatchCallback = std::function<void(index, const Match &)>;

    ~SubgraphIsomorphism() override = default;

    /**
     * Runs the search. If a callback throws or CTRL+C interrupts the search, run() stops all
     * workers and rethrows, and the results stay unavailable. Throws `std::runtime_error` if a
     * graph changed in a way that invalidates the input.
     */
    void run() override = 0;

    /// Number of workers @ref run() uses; 1 for the sequential algorithms.
    virtual count numberOfWorkers() const { return 1; }

    /**
     * Restricts matches to map every pattern node to a target node with the same label. The label
     * @ref none matches any label, and two empty vectors remove the labels. Call this before
     * @ref run().
     *
     * @param patternNodeLabels Labels of the pattern nodes, indexed by node id.
     * @param targetNodeLabels Labels of the target nodes, indexed by node id.
     */
    void setNodeLabels(const std::vector<index> &patternNodeLabels,
                       const std::vector<index> &targetNodeLabels);

    /**
     * Restricts matches to map every pattern edge to a target edge with the same label. The label
     * @ref none matches any label, and two empty vectors remove the labels. Both graphs need edge
     * ids, see `Graph::indexEdges()`. @ref run() throws for parallel edges with different labels.
     * Call this before @ref run().
     *
     * @param patternEdgeLabels Labels of the pattern edges, indexed by edge id.
     * @param targetEdgeLabels Labels of the target edges, indexed by edge id.
     */
    void setEdgeLabels(const std::vector<index> &patternEdgeLabels,
                       const std::vector<index> &targetEdgeLabels);

    /**
     * Passes every match to @a callback instead of storing it, so @ref getMatches() throws.
     * Replaces an earlier callback of either form. Call this before @ref run().
     */
    void setCallback(MatchCallback callback);

    /// Like @ref setCallback(MatchCallback), for a callback that may be called concurrently.
    void setCallback(ParallelMatchCallback callback);

    /**
     * Sets whether matches are stored. Pass false to only count them, so @ref getMatches() throws.
     * Call this before @ref run().
     */
    void setStoreMatches(bool storeMatches);

    /// Returns the matches of @ref run(). Throws if @ref run() has not finished, a callback was set
    /// or matches are not stored.
    const std::vector<Match> &getMatches() const;

    /**
     * Returns the number of matches found, stored or not. A parallel search with a callback may
     * deliver and count a few matches beyond `maxMatches`.
     */
    count numberOfMatches() const;

    /// Returns whether a match was found. Pass `maxMatches = 1` to stop at the first one.
    bool hasMatch() const;

protected:
    // A sequential run() calls prepareRun(), then reportMatch() per match until it returns false,
    // then finishRun(). A parallel run() calls prepareRun(), passes every match to
    // invokeCallback(), and calls finishRun(matches, found) after the join. Inside the parallel
    // region, workers may only poll Aux::SignalHandler::isRunning(). See ParallelRI.cpp.

    SubgraphIsomorphism(const Graph &pattern, const Graph &target, Semantics semantics,
                        count maxMatches);

    bool isNodeLabelled() const noexcept { return !patternNodeLabels.empty(); }

    /// An algorithm without edge-label support must throw from run() if this returns true.
    bool isEdgeLabelled() const noexcept { return !patternEdgeLabels.empty(); }

    bool hasCallback() const noexcept {
        return static_cast<bool>(callback) || static_cast<bool>(parallelCallback);
    }

    bool hasSerialCallback() const noexcept { return static_cast<bool>(callback); }

    bool storesMatches() const noexcept { return storeMatches && !hasCallback(); }

    /// Resets the results and revalidates the input. Call this at the start of run().
    void prepareRun();

    /// Records one match of a sequential search. Returns false once the cap is reached.
    bool reportMatch(const Match &match);

    /// How many more matches a sequential search may report, or @ref none without a cap.
    count remainingMatches() const noexcept {
        if (maxMatches == 0)
            return none;
        return matchCount >= maxMatches ? 0 : maxMatches - matchCount;
    }

    /**
     * Passes one match to the callback from any thread, but neither counts nor stores it. Returns
     * true if a callback received the match, in which case the caller must not store it.
     */
    bool invokeCallback(index tid, const Match &match);

    /// Marks a sequential run as finished.
    void finishRun();

    /**
     * Marks a parallel run as finished. Adopts the concatenated matches of all workers and, unless
     * a callback is set, trims them to the cap. Call this once, after the join.
     */
    void finishRun(std::vector<Match> &&matches, count found);

    const Graph *pattern;
    const Graph *target;

    std::vector<index> patternNodeLabels;
    std::vector<index> targetNodeLabels;

    std::vector<index> patternEdgeLabels;
    std::vector<index> targetEdgeLabels;

    Semantics semantics;
    /// 0 means no limit.
    count maxMatches;

private:
    void validateInput() const;

    void validateNodeLabels(const std::vector<index> &patternNodeLabels,
                            const std::vector<index> &targetNodeLabels) const;

    void validateEdgeLabels(const std::vector<index> &patternEdgeLabels,
                            const std::vector<index> &targetEdgeLabels) const;

    std::vector<Match> result;

    MatchCallback callback;
    ParallelMatchCallback parallelCallback;

    /// Serializes the calls of `callback`.
    std::mutex reportMutex;

    count matchCount;
    bool storeMatches;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_
