#ifndef NETWORKIT_ISOMORPHISM_VF2_HPP_
#define NETWORKIT_ISOMORPHISM_VF2_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {

/**
 * @ingroup isomorphism
 * Finds every occurrence of a pattern graph inside a target graph using the VF2 algorithm.
 *
 * VF2 extends a partial mapping one node pair at a time, depth first, and draws the candidate pairs
 * from the terminal sets of unmapped nodes adjacent to mapped ones. It supports both semantics,
 * directed and undirected graphs, node labels and edge labels. VF2 serves as the reference
 * implementation, but @ref RI is usually faster.
 *
 * The implementation is based on
 *
 * Cordella, L. P., Foggia, P., Sansone, C., & Vento, M. (2004).
 * A (Sub)Graph Isomorphism Algorithm for Matching Large Graphs.
 * IEEE Transactions on Pattern Analysis and Machine Intelligence, 26(10), 1367-1372.
 */
class VF2 final : public SubgraphIsomorphism {

public:
    /**
     * @param pattern The graph to look for. Must not contain self-loops.
     * @param target The graph to look in. Must agree with @a pattern on directedness.
     * @param semantics Whether matches must be induced. See @ref SubgraphIsomorphism::Semantics.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     */
    VF2(const Graph &pattern, const Graph &target, Semantics semantics = Semantics::INDUCED,
        count maxMatches = 0);

    void run() override;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_VF2_HPP_
