/*
 * Assortativity.hpp
 *
 *  Created on: Jun 13, 2015
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_CORRELATION_ASSORTATIVITY_HPP_
#define NETWORKIT_CORRELATION_ASSORTATIVITY_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup correlation
 *
 * Assortativity computes a coefficient that expresses the correlation of a
 * node attribute among connected pairs of nodes.
 *
 * For a partition of the node set (a nominal attribute), the coefficient is
 * computed from the mixing matrix e, where e_ij is the fraction of edge counts
 * from nodes in subset i to nodes in subset j. In undirected graphs, each edge
 * {u, v} with u != v is counted once in each direction and each self-loop is
 * counted once, and the mixing matrix is normalized by the total of these
 * counts (2m minus the number of self-loops). In directed graphs, each edge is
 * counted once and the matrix is normalized by m.
 */
class Assortativity final : public Algorithm {

public:
    /**
     * Initialize Assortativity with a graph @a G and an array of numerical
     * node values.
     *
     * @param G The graph.
     * @param attribute		numerical node value array
     */
    Assortativity(const Graph &G, const std::vector<double> &attribute);

    /**
     * Initialize Assortativity with a graph @a G and a partition of the node set
     *
     * @param G The graph.
     * @param partition		partition of the node set
     */
    Assortativity(const Graph &G, const Partition &partition);

    /**
     * Runs the algorithm. The algorithm is not parallel.
     */
    void run() override;

    /**
     * Return the assortativity coefficient.
     */
    double getCoefficient() const;

private:
    const Graph *G;
    const std::vector<double> emptyVector;
    const Partition emptyPartition;
    const std::vector<double> *attribute;
    const Partition *partition;
    bool nominal; // whether we calculate assortativity for a nominal or ordinal attribute
    double coefficient;
};

} /* namespace NetworKit */
#endif // NETWORKIT_CORRELATION_ASSORTATIVITY_HPP_
