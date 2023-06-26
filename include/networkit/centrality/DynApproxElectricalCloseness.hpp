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
#include <networkit/centrality/Centrality.hpp>
#include <networkit/components/BiconnectedComponents.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class DynApproxElectricalCloseness final : public Centrality, public DynAlgorithm {

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
     * This dynamic algorithm supports addition of arbitrary edges and deletion of edges which are
     * not in the bfs tree sourced from the pivot node. The algorithm depends on G having a single
     * connected component - edge deletions which disconnect the graph are not supported.
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

    void updateBatch([[maybe_unused]] const std::vector<GraphEvent> &batch) override {
        throw std::runtime_error("Error: batch updates are not supported.");
    }

    /**
     * Return an epsilon-approximation of the diagonal of the laplacian's pseudoinverse.
     *
     * @return Approximation of the diagonal of the laplacian's pseudoinverse.
     */
    const std::vector<double> &getDiagonal() const {
        assureFinished();
        return diagonal;
    }

    inline NetworKit::count getUstCount() { return numberOfUSTs; }

private:
    const double epsilon, delta, kappa;
    double tol;
    count numberOfUSTs;
    node root = 0;
    const node pivot;
    count rootEcc;
    NetworKit::CSRMatrix laplacian;

    // #of BFSs used to estimate a vertex with low eccentricity.
    static constexpr uint32_t sweeps = 10;

    enum class NodeStatus : unsigned char {
        NOT_IN_COMPONENT,
        IN_SPANNING_TREE,
        NOT_VISITED,
        VISITED
    };

    struct Tree {
        std::vector<node> parent;
        std::vector<node> sibling;
        std::vector<node> child;

        std::vector<count> tVisit;
        std::vector<count> tFinish;
        bool timesComputed = false;
    };

    // Used to mark the status of each node, one vector per thread
    std::vector<std::vector<NodeStatus>> statusGlobal;

    BiconnectedComponents bcc;

    // Nodes in each biconnected components sorted by their degree.
    std::vector<std::vector<node>> sequences;

    // Shortest path tree
    Tree bfsTree;

    // Index of the parent component of the current component (after the topological order has been
    // determined)
    std::vector<index> biParent;
    // Node within the bionnected component that points to the node in the parent component
    std::vector<node> biAnchor;

    // Topological order of the biconencted components
    std::vector<index> topOrder;

    // approx effective resistance, summed contribution per tree, one per thread.
    std::vector<std::vector<double>> approxEffResistanceGlobal;

    // Pseudo-inverse diagonal
    std::vector<double> diagonal;

    // Resistance of node with root
    std::vector<double> resistanceToRoot;

    // lpinv column of root
    Vector rootCol;

    // Least common ancestor vectors, per thread
    std::vector<std::vector<node>> lcaGlobal;

    // Random number generators
    std::vector<std::mt19937_64> generators;
    std::vector<std::uniform_int_distribution<index>> degDist;

    // Nodes sequences: Wilson's algorithm runs faster if we start the random walks following a
    // specific sequence of nodes. In this function we compute those sequences.
    void computeNodeSequence();

    void setDFSTimes(Tree &tree);

    void computeBFSTree();
    void sampleUST(Tree &result);
    void sampleUSTWithEdge(Tree &result, node a, node b);

    void aggregateUST(Tree &tree, double weight);

    void edgeAdded(node a, node b);

    void edgeRemoved(node a, node b);

    node approxMinEccNode();

    count computeNumberOfUSTs() const;

#ifdef NETWORKIT_SANITY_CHECKS
    // Debugging methods
    void checkBFSTree() const;
    void checkUST(const Tree &tree) const;
    void checkTwoNodesSequence(const std::vector<node> &sequence, std::vector<node> &parent) const;
    void checkTimeStamps(const Tree &tree) const;
#endif // NETWORKIT_SANITY_CHECKS
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_DYN_APPROX_ELECTRICAL_CLOSENESS_HPP_
