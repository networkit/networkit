#ifndef NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_
#define NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_

#include <functional>
#include <mutex>
#include <span>
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
 * edges. The algorithm holds both graphs by reference, so they must outlive it. Since the graphs
 * may change after construction, @ref run() validates them and the labels each time.
 */
class SubgraphIsomorphism : public Algorithm {

public:
    // The underlying type stays int. Python passes a C int, and a narrower type would wrap an
    // invalid value such as 256 to a valid one before the constructor could reject it.
    enum class Semantics {
        /// Pattern edges map to target edges and pattern non-edges map to target non-edges.
        INDUCED,
        /// Pattern edges map to target edges. Additional target edges are allowed.
        MONOMORPHISM,
    };

    /// One match: `match[u]` is the image of pattern node u, or @ref none if u is not a node.
    using Match = std::vector<node>;

    /**
     * Receives one match at a time, never concurrently. The match is indexed like a @ref Match.
     * It refers to an internal buffer and is valid only during the call. A parallel algorithm
     * serializes the calls with a lock, which a @ref ParallelMatchCallback avoids.
     */
    using MatchCallback = std::function<void(std::span<const node>)>;

    /**
     * Receives one match at a time together with the id of the worker that found it. The id is
     * less than `Aux::getMaxNumberOfThreads()` at the start of @ref run(). The callback may be
     * called concurrently and must be thread-safe. The sequential algorithms pass worker id 0.
     */
    using ParallelMatchCallback = std::function<void(index, std::span<const node>)>;

    ~SubgraphIsomorphism() override = default;

    /**
     * Runs the search. If a callback throws or CTRL+C interrupts the search, run() stops all
     * workers and rethrows, and the results stay unavailable.
     *
     * @throws std::runtime_error If the graphs break the requirements above, or if the labels do
     * not fit the graphs.
     */
    void run() override = 0;

    /**
     * Restricts matches to map every pattern node to a target node with the same label. The label
     * @ref none matches any label, and two empty vectors remove the labels. The labels are copied.
     * Call this before @ref run(), which throws if a vector is shorter than the node id bound of
     * its graph.
     *
     * @param patternNodeLabels Labels of the pattern nodes, indexed by node id.
     * @param targetNodeLabels Labels of the target nodes, indexed by node id.
     */
    void setNodeLabels(std::span<const index> patternNodeLabels,
                       std::span<const index> targetNodeLabels);

    /**
     * Restricts matches to map every pattern edge to a target edge with the same label. The label
     * @ref none matches any label, and two empty vectors remove the labels. The labels are copied.
     * Call this before @ref run(), which throws if a graph has no edge ids (see
     * `Graph::indexEdges()`), if a vector is shorter than the edge id bound of its graph, or for
     * parallel edges with different labels.
     *
     * @param patternEdgeLabels Labels of the pattern edges, indexed by edge id.
     * @param targetEdgeLabels Labels of the target edges, indexed by edge id.
     */
    void setEdgeLabels(std::span<const index> patternEdgeLabels,
                       std::span<const index> targetEdgeLabels);

    /**
     * Passes every match to @a callback instead of storing it, so @ref getMatches() throws.
     * Replaces an earlier callback of either form. Call this before @ref run().
     */
    void setCallback(MatchCallback callback);

    /// Like @ref setCallback(), for a callback that may be called concurrently.
    void setParallelCallback(ParallelMatchCallback callback);

    /**
     * Sets whether matches are stored. Pass false to only count them, so @ref getMatches() throws.
     * Call this before @ref run().
     */
    void setStoreMatches(bool storeMatches);

    /**
     * Returns the matches of the last @ref run(). Throws if @ref run() has not finished, or if that
     * run did not store the matches because a callback was set or matches were only counted.
     */
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

    /// @throws std::invalid_argument If @a semantics is not a known value.
    SubgraphIsomorphism(const Graph &pattern, const Graph &target, Semantics semantics,
                        count maxMatches);

    bool hasCallback() const noexcept {
        return static_cast<bool>(callback) || static_cast<bool>(parallelCallback);
    }

    /// Whether the current run stores its matches. prepareRun() fixes the value from the settings.
    bool storesMatches() const noexcept { return matchesStored; }

    /// Resets the results and validates the input. Call this at the start of run().
    void prepareRun();

    /// Records one match of a sequential search. Returns false once the cap is reached.
    bool reportMatch(std::span<const node> match);

    /**
     * Passes one match to the callback from any thread, but neither counts nor stores it. Returns
     * true if a callback received the match, in which case the caller must not store it.
     */
    bool invokeCallback(index tid, std::span<const node> match);

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
    std::vector<Match> result;

    MatchCallback callback;
    ParallelMatchCallback parallelCallback;

    /// Serializes the calls of `callback`.
    std::mutex reportMutex;

    count matchCount;
    bool storeMatches;
    /// Whether the last run() stores its matches. The settings may change after that run, so
    /// getMatches() checks this flag rather than the settings.
    bool matchesStored;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_SUBGRAPH_ISOMORPHISM_HPP_
