#ifndef NETWORKIT_ISOMORPHISM_RI_HPP_
#define NETWORKIT_ISOMORPHISM_RI_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {

/**
 * @ingroup isomorphism
 * Finds every occurrence of a pattern graph inside a target graph using the RI algorithm.
 *
 * RI fixes the order of the pattern nodes once and then backtracks along this order, without the
 * terminal sets of @ref VF2. It supports both semantics, directed and undirected graphs, node
 * labels and edge labels.
 *
 * The variant @ref Variant::RI_DS first computes a domain of candidate target nodes for every
 * pattern node. It pays off for disconnected patterns and selective node labels, but otherwise
 * usually costs more than it saves.
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
    // The underlying type stays int. Python passes a C int, and a narrower type would wrap an
    // invalid value such as 256 to a valid one before the constructor could reject it.
    enum class Variant {
        /// Plain RI.
        RI,
        /// RI-DS-SI-FC: candidate domains with forward checking.
        RI_DS
    };

    /**
     * @param pattern The graph to look for. Must not contain self-loops.
     * @param target The graph to look in. Must agree with @a pattern on directedness.
     * @param variant Plain RI or RI-DS. See @ref Variant.
     * @param semantics Whether matches must be induced. See @ref SubgraphIsomorphism::Semantics.
     * @param maxMatches Stop after this many matches; 0 means no limit.
     * @param domainFilterThreshold Only for @ref Variant::RI_DS, which builds the domains from
     * degrees and labels and then shrinks them. The search filters the candidates of a pattern
     * node through its domain if the shrinking removed at least this fraction of the domain or left
     * one target node. Below the threshold, the filter removes too few candidates to pay for its
     * lookups. 0 always filters, and 1 filters only through single-node domains. The value changes
     * the running time, not the matches. The default comes from caidaRouterLevel with a triangle
     * pattern and node labels, where filtering cost time below about 0.8.
     *
     * @throws std::invalid_argument If @a variant or @a semantics is not a known value, or if
     * @a domainFilterThreshold is not between 0 and 1. @ref run() checks the graphs, since they may
     * change after construction.
     */
    RI(const Graph &pattern, const Graph &target, Variant variant = Variant::RI,
       Semantics semantics = Semantics::INDUCED, count maxMatches = 0,
       double domainFilterThreshold = 0.8);

    void run() override;

private:
    Variant variant;
    double domainFilterThreshold;
};

} // namespace NetworKit

#endif // NETWORKIT_ISOMORPHISM_RI_HPP_
