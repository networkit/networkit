/*
 * ApproxElectricalCloseness.hpp
 *
 *  Created on: 17.10.2019
 *     Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *              Alexander van der Grinten <avdgrinten@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_APPROX_ELECTRICAL_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_APPROX_ELECTRICAL_CLOSENESS_HPP_

#include <cmath>
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/components/BiconnectedComponents.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class ApproxElectricalCloseness : public Centrality {

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
     * @param G The input graph.
     * @param epsilon Maximum absolute error of the elements in the diagonal.
     * @param kappa Balances the tolerance of the solver for the linear system and the number of
     * USTs to be sampled.
     */
    ApproxElectricalCloseness(const Graph &G, double epsilon = 0.1, double kappa = 0.3);

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
     * @param G The input graph.
     * @param epsilon Maximum absolute error of the elements in the diagonal.
     * @param kappa Balances the tolerance of the solver for the linear system and the number of
     * USTs to be sampled.
     * @param delta approximation probability.
     */
    ApproxElectricalCloseness(const Graph &G, double epsilon, double kappa, double delta);

    ~ApproxElectricalCloseness() override = default;

    /**
     * Run the algorithm.
     */
    void run() override;

    /**
     * Return an epsilon-approximation of the diagonal of the laplacian's pseudoinverse.
     *
     * @return Approximation of the diagonal of the laplacian's pseudoinverse.
     */
    const std::vector<double> &getDiagonal() const {
        assureFinished();
        return diagonal;
    }

    /**
     * Compute and return the nearly-exact values of the diagonal of the laplacian's pseudoinverse.
     * The values are computed by solving Lx = e_u - 1 / n for every vertex u of the graph with a
     * LAMG solver.
     *
     * @param tol Tolerance for the LAMG solver.
     *
     * @return Nearly-exact values of the diagonal of the laplacian's pseudoinverse.
     */
    std::vector<double> computeExactDiagonal(double tol = 1e-9) const;

protected:
    const double epsilon, delta, kappa;
    double tol;
    node root = 0;
    count rootEcc;
    count numberOfUSTs;

    // These members are artifacts from running the algorithm. In the static algorithm, these
    // members are cleared after running the algorithm to release memory; in the dynamic algorithm,
    // they are kept to allow updates.

    // laplacian matrix of the graph.
    CSRMatrix laplacian;
    // lpinv column of root
    Vector rootCol;
    // Resistance of node with root
    std::vector<double> resistanceToRoot;

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

        Tree(count n) : parent(n, none), sibling(n), child(n), tVisit(n), tFinish(n) {}
    };

    std::vector<Tree> bfsTrees;

    // Used to mark the status of each node, one vector per thread
    std::vector<std::vector<NodeStatus>> statusGlobal;

    BiconnectedComponents bcc;

    // Nodes in each biconnected components sorted by their degree.
    std::vector<std::vector<node>> sequences;

    // Index of the parent component of the current component (after the topological order has been
    // determined)
    std::vector<index> biParent;
    // Node within the biconnected component that points to the node in the parent component
    std::vector<node> biAnchor;

    // Topological order of the biconnected components
    std::vector<index> topOrder;

    // normalized approx effective resistance, one per thread.
    std::vector<std::vector<double>> approxEffResistanceGlobal;

    // Pseudo-inverse diagonal
    std::vector<double> diagonal;

    // Random number generators
    std::vector<std::mt19937_64> generators;
    std::vector<std::uniform_int_distribution<index>> degDist;

    // Parent pointers of the bfs tree
    std::vector<node> bfsParent;

    // Nodes sequences: Wilson's algorithm runs faster if we start the random walks following a
    // specific sequence of nodes. In this function we compute those sequences.
    void computeNodeSequence(node pivot = none);

    void computeBFSTree();
    void sampleUST();
    void dfsUST();
    void aggregateUST(double weight = 1.0);

    node approxMinEccNode();

    void aggregateResults();
    void solveColumnAndSampleUSTs();

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

#endif // NETWORKIT_CENTRALITY_APPROX_ELECTRICAL_CLOSENESS_HPP_
