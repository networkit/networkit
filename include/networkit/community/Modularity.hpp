/*
 * Modularity.hpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COMMUNITY_MODULARITY_HPP_
#define NETWORKIT_COMMUNITY_MODULARITY_HPP_

#include <networkit/community/QualityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Modularity is a quality index for community detection. It assigns a quality value in [-0.5, 1.0]
 * to a partition of a graph which is higher for more modular networks and partitions which better
 * capture the modular structure.
 *
 * Modularity is defined as
 *
 * $$mod(\\zeta) := cov(\\zeta)
 *     - \\sum_{C \\in \\zeta} \\frac{vol^{out}(C) \\cdot vol^{in}(C)}{\\omega(E)^2},$$
 *
 * where $\\omega(E)$ is the total edge weight. For directed graphs (Leicht and Newman, 2008),
 * $vol^{out}(C)$ and $vol^{in}(C)$ are the out- and in-volumes of $C$. For undirected graphs
 * both are half the volume of $C$.
 */
class Modularity final : public QualityMeasure {

private:
    double gTotalEdgeWeight;

public:
    /** Default constructor */
    Modularity();

    /**
     * Returns the Modularity of the given clustering with respect to the graph @a G.
     *
     * @param zeta The clustering.
     * @param G The graph.
     * @return The modularity.
     */
    double getQuality(const Partition &zeta, const Graph &G) override;

    /**
     * @param totalEdgeWeight Sum of all edge weights in @a G. If specified, it does not
     *        have to be computed.
     */
    void setTotalEdgeWeight(double totalEdgeWeight);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_MODULARITY_HPP_
