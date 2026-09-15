#ifndef NETWORKIT_ISOMORPHISM_RI_HPP_
#define NETWORKIT_ISOMORPHISM_RI_HPP_

#include <cstdint>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {

/**
 * @ingroup isomorphism
 * Finds every occurrence of a pattern graph inside a target graph using the RI algorithm.
 *
 * See @ref SubgraphIsomorphism for the definition of a match and how to use the class.
 *
 * RI computes a fixed order of the pattern nodes once and then runs a plain backtracking search
 * along this order, without the terminal sets of @ref VF2. The order starts at a node of maximum
 * degree and repeatedly appends the unordered node with the most edges into the ordered nodes.
 * Ties are broken first by the number of ordered nodes reachable in two hops through an unordered
 * node, then by the number of unordered neighbours that are adjacent to no ordered node, and
 * finally by the smallest node id. During the search, the candidates for a pattern node are the
 * neighbours of the target node that an earlier adjacent pattern node was mapped to. A pattern
 * node without an earlier neighbour takes every target node as a candidate.
 *
 * The variant @ref Variant::RI_DS computes a domain of candidate target nodes for every pattern
 * node before the search, from degrees, node labels and one refinement pass over the pattern
 * edges. Pattern nodes whose domain holds a single target node are ordered first, and forward
 * checking removes that target node from all other domains. Remaining ties in the order are broken
 * by domain size before node id. RI-DS pays off for disconnected patterns and for selective node
 * labels. Otherwise, computing the domains usually costs more than it saves, so @ref Variant::RI
 * is the default.
 *
 * RI supports both semantics, directed and undirected graphs, node labels and edge labels. The
 * search can be interrupted with CTRL+C.
 *
 * The implementation is based on
 *
 * Bonnici, V., Giugno, R., Pulvirenti, A., Shasha, D., & Ferro, A. (2013).
 * A subgraph isomorphism algorithm and its application to biochemical data.
 * BMC Bioinformatics, 14(Suppl 7), S13.
 *
 * and, for the RI-DS variant, which the paper calls RI-DS-SI-FC, on
 *
 * Kimmig, R., Meyerhenke, H., & Strash, D. (2017).
 * Shared Memory Parallel Subgraph Enumeration.
 * IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW).
 */
class RI final : public SubgraphIsomorphism {

public:
    /**
     * Which variant of RI to run.
     */
    enum class Variant : uint8_t {
        /// Plain RI. The default.
        RI,
        /// RI-DS-SI-FC: candidate domains with forward checking, and an order that puts pattern
        /// nodes with a single candidate first and breaks ties by domain size.
        RI_DS
    };

    /**
     * @param pattern The graph to look for. Must not contain self-loops.
     * @param target The graph to look in. Must agree with @a pattern on directedness.
     * @param variant Plain RI or RI-DS. See @ref Variant.
     * @param semantics Whether matches must be induced. See @ref SubgraphIsomorphism::Semantics.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     */
    RI(const Graph &pattern, const Graph &target, Variant variant = Variant::RI,
       Semantics semantics = Semantics::INDUCED, count maxMatches = 0);

    /**
     * Runs the search. Query the result with @ref getMatches(), @ref numberOfMatches() or
     * @ref hasMatch().
     */
    void run() override;

private:
    Variant variant;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_RI_HPP_
