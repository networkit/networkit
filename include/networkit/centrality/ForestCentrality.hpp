/*
 * ForestCentrality.hpp
 *
 *  Created on: 12.02.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_FOREST_CENTRALITY_HPP_
#define NETWORKIT_CENTRALITY_FOREST_CENTRALITY_HPP_

#include <cmath>
#include <functional>
#include <omp.h>
#include <random>
#include <vector>

#include <networkit/algebraic/Vector.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

// NOTE: Does not support multiple calls of run()
class ForestCentrality final : public Centrality {

public:
    /**
     * Approximates the forest closeness centrality of all the vertices of a graph by approximating
     * the diagonal of the forest matrix of @a G. Every element of the diagonal is guaranteed to
     * have a maximum absolute error of @a epsilon. Based on "New Approximation Algorithms for
     * Forest Closeness Centrality - for Individual Vertices and Vertex Groups", van der Grinten et
     * al, SDM 2021.
     * The algorithm runs in two steps: (i) solving a linear system and (ii) sampling uniform
     * spanning trees (USTs). The parameter @a kappa balances the tolerance of the linear sytem
     * solver and the number of USTs to be sampled. A high value of @a kappa raises the tolerance
     * (solver converges faster) but more USTs need to be sampled, vice versa for a low value of @a
     * kappa. Note: the algorithm requires an augmented graph as input (see the reference paper for
     * details). An augmented graphs can be generated with GraphTools::createAugmentedGraph.
     *
     * @param G The input graph. Must be an augmented graph (an augmented graph can be crated with
     * GraphTools::createAugmentedGraph.
     * @param root Root node of the augmented graph.
     * @param epsilon Maximum absolute error of the elements in the diagonal.
     * @param kappa Balances the tolerance of the linear system solver and the number of USTs to be
     * sampled.
     */
    explicit ForestCentrality(const Graph &G, node root, double epsilon = 0.1, double kappa = 0.3);

    /**
     * Run the algorithm.
     */
    void run() override;

    /**
     * Return an epsilon-approximation of the diagonal of the forest matrix.
     *
     * @return Approximation of the diagonal of the forest matrix.
     */
    const std::vector<double> &getDiagonal() const {
        assureFinished();
        return diagonal;
    }

    /**
     * Return the number of sampled USTs.
     *
     * @return Number of sampled USTs.
     */
    count getNumberOfSamples() const noexcept { return numberOfUSTs; }

private:
    node root;
    const double epsilon, kappa, volG;
    const count numberOfUSTs;

    // Used to mark the status of each node, one vector per thread
    std::vector<std::vector<uint8_t>> statusGlobal;

    // Pointers to the parent of the UST, one vector per thread
    std::vector<std::vector<node>> parentsGlobal;

    // Nodes sorted by decreasing degree
    std::vector<node> decDegree;

    // Non-normalized approx effective resistance, one per thread.
    std::vector<std::vector<count>> approxEffResistanceGlobal;

    // Distributions for selecting random neighbors
    std::vector<std::uniform_int_distribution<count>> uniformDistr;

    // Solution of the linear system
    Vector linearSysSol;

    // Diagonal elements
    std::vector<double> diagonal;

    count computeNumberOfUSTs() const {
        return count{4}
               * static_cast<count>(std::ceil(
                   std::log(2.0 * static_cast<double>(G.numberOfEdges()) * G.numberOfNodes())
                   / (2.0 * epsilon * epsilon)));
    }

    void sampleUSTs();
    void solveLinearSystem();
    void computeDiagonal();
    void computeScores();
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_FOREST_CENTRALITY_HPP_
