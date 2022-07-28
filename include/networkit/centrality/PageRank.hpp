/*
 * PageRank.h
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_CENTRALITY_PAGE_RANK_HPP_
#define NETWORKIT_CENTRALITY_PAGE_RANK_HPP_

#include <atomic>
#include <limits>
#include <memory>

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Compute PageRank as node centrality measure. In the default mode this computation is in line
 * with the original paper "The PageRank citation ranking: Bringing order to the web." by L. Brin et
 * al (1999). In later publications ("PageRank revisited." by M. Brinkmeyer et al. (2005) amongst
 * others) sink-node handling was added for directed graphs in order to comply with the theoretical
 * assumptions by the underlying Markov chain model. This can be activated by setting the matching
 * parameter to true. By default this is disabled, since it is an addition to the original
 * definition.
 *
 * Page-Rank values can also be normalized by post-processed according to "Comparing Apples and
 * Oranges: Normalized PageRank for Evolving Graphs" by Berberich et al. (2007). This decouples
 * the PageRank values from the size of the input graph. To enable this, set the matching parameter
 * to true. Note that sink-node handling is automatically activated if normalization is used.
 *
 * NOTE: There is an inconsistency in the definition in Newman's book (Ch. 7) regarding
 * directed graphs; we follow the verbal description, which requires to sum over the incoming
 * edges (as opposed to outgoing ones).
 */
class PageRank final : public Centrality {

public:
    enum Norm {
        L1_NORM,
        L2_NORM,
        L1Norm = L1_NORM, // this + following added for backwards compatibility
        L2Norm = L2_NORM
    };

    enum SinkHandling { NO_SINK_HANDLING, DISTRIBUTE_SINKS };

    /**
     * Constructs the PageRank class for the Graph @a G
     *
     * @param[in] G Graph to be processed.
     * @param[in] damp Damping factor of the PageRank algorithm.
     * @param[in] tol Error tolerance for PageRank iteration.
     * @param[in] distributeSinks Set to distribute PageRank values for sink nodes. (default =
     * false)
     * @param[in] normalized Set to normalize Page-Rank values in order to decouple it from the
     * network size. (default = false)
     */
    PageRank(const Graph &G, double damp = 0.85, double tol = 1e-8, bool normalized = false,
             SinkHandling distributeSinks = SinkHandling::NO_SINK_HANDLING);

    /**
     * Computes page rank on the graph passed in constructor.
     */
    void run() override;

    /**
     * Returns the maximum PageRank score
     */
    double maximum() override;

    /**
     * Return the number of iterations performed by the algorithm.
     *
     * @return Number of iterations performed by the algorithm.
     */
    count numberOfIterations() const {
        assureFinished();
        return iterations;
    }

    // Maximum number of iterations allowed
    count maxIterations = std::numeric_limits<count>::max();

    // Norm used as stopping criterion
    Norm norm = Norm::L2_NORM;

private:
    double damp;
    double tol;
    count iterations;
    bool normalized;
    SinkHandling distributeSinks;
    std::atomic<double> max;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_PAGE_RANK_HPP_
