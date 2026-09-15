#ifndef NETWORKIT_ISOMORPHISM_VF3_HPP_
#define NETWORKIT_ISOMORPHISM_VF3_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {

/**
 * @ingroup isomorphism
 * Finds every occurrence of a pattern graph inside a target graph using the VF3 algorithm.
 *
 * @warning The search is not implemented yet, so @ref run() throws. Use @ref RI or @ref VF2
 * instead.
 *
 * See @ref SubgraphIsomorphism for the definition of a match and how to use the class.
 *
 * VF3 extends @ref VF2 in three ways. Nodes with the same label form a class, and the search
 * reasons about classes instead of single nodes. The order in which pattern nodes are mapped is
 * computed once before the search and prefers rare classes and nodes with many edges into the
 * ordered nodes. The feasibility counts are precomputed per class, so the pruning rules become
 * table lookups. VF3 is designed for large and dense targets, especially with node labels.
 *
 * VF3 supports both semantics, directed and undirected graphs, and node labels. It does not
 * support edge labels.
 *
 * VF3 is described in
 *
 * Carletti, V., Foggia, P., Saggese, A., & Vento, M. (2018).
 * Challenging the Time Complexity of Exact Subgraph Isomorphism
 * for Huge and Dense Graphs with VF3.
 * IEEE Transactions on Pattern Analysis and Machine Intelligence, 40(4), 804-818.
 */
class VF3 final : public SubgraphIsomorphism {

public:
    /**
     * @param pattern The graph to look for. Must not contain self-loops.
     * @param target The graph to look in. Must agree with @a pattern on directedness.
     * @param semantics Whether matches must be induced. See @ref SubgraphIsomorphism::Semantics.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     */
    VF3(const Graph &pattern, const Graph &target, Semantics semantics = Semantics::INDUCED,
        count maxMatches = 0);

    /**
     * Runs the search.
     *
     * @warning Not implemented. Throws `std::runtime_error` if edge labels were set with
     * @ref setEdgeLabels() or the input is invalid, and `std::logic_error` otherwise.
     */
    void run() override;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_VF3_HPP_
