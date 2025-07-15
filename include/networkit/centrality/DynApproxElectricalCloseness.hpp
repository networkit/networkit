/*
 * DynApproxElectricalCloseness.hpp
 *
 *  Created on: 19.04.2023
 *     Authors: Matthias Görg <goergmat@informatik.hu-berlin.de>
 *              Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_DYN_APPROX_ELECTRICAL_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_DYN_APPROX_ELECTRICAL_CLOSENESS_HPP_

#include <cmath>
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/components/BiconnectedComponents.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class DynApproxElectricalCloseness final : public ApproxElectricalCloseness, public DynAlgorithm {

public:
    /**
     * Approximates the electrical closeness of all the vertices of the graph by approximating the
     * diagonal of the laplacian's pseudoinverse of @a G. Every element of the diagonal has
     * a maximum absolute error of @a epsilon with probability (1-(1/n)). Based on "Approximation of
     * the Diagonal of a Laplacian’s Pseudoinverse for Complex Network Analysis", Angriman et al.,
     * ESA 2020. The algorithm does two steps: solves a linear system and samples uniform spanning
     * trees (USTs). The parameter @a kappa balances the tolerance of solver for the linear system
     * and the number of USTs to be sampled. A high value of @a kappa raises the tolerance (solver
     * converges faster) but more USTs need to be sampled, vice versa for a low value of @a kappa.
     *
     * @note This dynamic algorithm supports addition of arbitrary edges and deletion of edges which
     * are not in the bfs tree sourced from the pivot node. The algorithm depends on G having a
     * single connected component - edge deletions which disconnect the graph are not supported.
     * @note Batch updates are not supported.
     *
     * @param G The input graph.
     * @param epsilon Maximum absolute error of the elements in the diagonal.
     * @param kappa Balances the tolerance of the solver for the linear system and the number of
     * @param pivot The pivot node. If none, a node with small approximate
     * eccentricity will be chosen. USTs to be sampled.
     */
    DynApproxElectricalCloseness(const Graph &G, double epsilon = 0.1, double kappa = 0.3,
                                 node pivot = none);

    /**
     * Constructor to set @a delta (the approximation probability) manually.
     */
    DynApproxElectricalCloseness(const Graph &G, double epsilon, double kappa, node pivot,
                                 double delta);

    /** copy constructor that copies internal data structures, but still references the same graph
     * as @a other */
    DynApproxElectricalCloseness(const DynApproxElectricalCloseness &other) = default;

    ~DynApproxElectricalCloseness() override = default;

    /**
     * Run the algorithm.
     */
    void run() override;

    /**
     * update - to be called after the graph has been updated.
     * Supports arbitrary edge additions.
     * Supports edge deletions for edges which are not on the bfs tree of the pivot node. Removing
     * the edge may not disconnect the graph.
     */
    void update(GraphEvent e) override;

    /**
     * @note batch updates are not supported.
     */
    void updateBatch(const std::vector<GraphEvent> &) override {
        throw std::runtime_error("Error: batch updates are not supported.");
    }

    inline NetworKit::count getUstCount() { return numberOfUSTs; }

private:
    const node pivot;

    // Least common ancestor vectors, per thread
    std::vector<std::vector<node>> lcaGlobal;

    void sampleUSTWithEdge(node a, node b);

    void edgeAdded(node a, node b);
    void edgeRemoved(node a, node b);

    std::tuple<double, Vector, Vector> computeColumnsExact(node a, node b);
    void solveRootAndResample(node a, node b, double w);
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_DYN_APPROX_ELECTRICAL_CLOSENESS_HPP_
