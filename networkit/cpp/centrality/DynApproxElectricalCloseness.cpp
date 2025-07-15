/*
 * DynApproxElectricalCloseness.cpp
 *
 *  Created on: 19.04.2023
 *     Authors: Matthias Görg <goergmat@informatik.hu-berlin.de>
 *              Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */

#include <algorithm>
#include <numeric>
#include <omp.h>
#include <queue>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/DynApproxElectricalCloseness.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>
#include <networkit/numerics/Preconditioner/DiagonalPreconditioner.hpp>

#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

DynApproxElectricalCloseness::DynApproxElectricalCloseness(const Graph &G, double epsilon,
                                                           double kappa, node pivot, double delta)
    : ApproxElectricalCloseness(G, epsilon, kappa, delta), pivot(pivot) {

    const count n = G.upperNodeIdBound(), threads = omp_get_max_threads();
    lcaGlobal.resize(threads, std::vector<node>(n, none));
}

DynApproxElectricalCloseness::DynApproxElectricalCloseness(const Graph &G, double epsilon,
                                                           double kappa, node pivot)
    : DynApproxElectricalCloseness(G, epsilon, kappa, pivot,
                                   1.0 / static_cast<double>(G.numberOfNodes())) {}

void DynApproxElectricalCloseness::sampleUSTWithEdge(node a, node b) {
    auto n = G.numberOfNodes();
    assert(a < n && b < n);

    auto &tree = bfsTrees[omp_get_thread_num()];

    auto &parent = tree.parent;
    auto &status = statusGlobal[omp_get_thread_num()];
    auto &childPtr = tree.child;
    auto &siblingPtr = tree.sibling;

    std::fill(status.begin(), status.end(), NodeStatus::NOT_IN_COMPONENT);
    std::fill(parent.begin(), parent.end(), none);
    std::fill(childPtr.begin(), childPtr.end(), none);
    std::fill(siblingPtr.begin(), siblingPtr.end(), none);

    auto &generator = generators[omp_get_thread_num()];

    // Initialize UST with only (a, b) in the tree.
    parent[a] = b;
    status[a] = NodeStatus::IN_SPANNING_TREE;
    status[b] = NodeStatus::IN_SPANNING_TREE;
    unsigned int nodesInSpanningTree = 2;

    // The tree is generated using Wilson's algorithm, rooted in b.
    // Afterwards reroot it to this->root.
    for (const auto componentIndex : topOrder) {
        const auto &sequence = sequences[componentIndex];
        for (const auto walkStart : sequence) {
            if (status[walkStart] == NodeStatus::IN_SPANNING_TREE) {
                continue;
            }

            // Random Walk
            node currentNode = walkStart;
            do {
                node randomNeighbor =
                    G.getIthNeighbor(currentNode, degDist[currentNode](generator));

                assert(randomNeighbor != none);
                parent[currentNode] = randomNeighbor;
                currentNode = randomNeighbor;

            } while (status[currentNode] != NodeStatus::IN_SPANNING_TREE);

            // Last node encountered in the random walk (it is in the spanning tree);
            const node walkEnd = currentNode;
            assert(status[walkEnd] == NodeStatus::IN_SPANNING_TREE);
            // Add the random walk to the spanning tree;
            for (currentNode = walkStart; currentNode != walkEnd;
                 currentNode = parent[currentNode]) {

                status[currentNode] = NodeStatus::IN_SPANNING_TREE;
                ++nodesInSpanningTree;
            }
        }
    }

    assert(nodesInSpanningTree == n);

    // Switch root from b to `root`.
    node currentNode = root;
    node next = parent[root];
    while (next != none) {
        node tmp = parent[next];
        parent[next] = currentNode;
        currentNode = next;
        next = tmp;
    }
    parent[root] = none;

    // Set child and sibling pointers
    count visitedNodes = 0;
    for (node u : G.nodeRange()) {
        while (status[u] == NodeStatus::IN_SPANNING_TREE) {
            status[u] = NodeStatus::NOT_VISITED;
            ++visitedNodes;
            node parentU = parent[u];
            if (parentU != none) {
                assert(siblingPtr[u] == none);
                if (childPtr[parentU] != none)
                    siblingPtr[u] = childPtr[parentU];
                childPtr[parentU] = u;
                u = parentU;
            } else
                break;
        }
        if (visitedNodes == G.numberOfNodes())
            break;
    }

#ifdef NETWORKIT_SANITY_CHECKS
    checkUST(tree);
#endif
}

void DynApproxElectricalCloseness::edgeAdded(node a, node b) {
    assert(G.hasEdge(a, b));
    assert(hasRun);

    const count n = G.numberOfNodes();
    const double n_double = static_cast<double>(n);

    // Compute lpinv columns for a and b exactly.

    laplacian.setValue(a, a, laplacian(a, a) + 1.);
    laplacian.setValue(b, b, laplacian(b, b) + 1.);
    laplacian.setValue(a, b, laplacian(a, b) - 1.);
    laplacian.setValue(b, a, laplacian(b, a) - 1.);

    const auto [w, lpinvColA, lpinvColB] = computeColumnsExact(a, b);

    solveRootAndResample(a, b, w);

    // Aggregate results
    Vector currentRoundResistanceApprox(n);

    G.parallelForNodes([&](const node u) {
        // Accumulate all results on first thread vector
        for (count i = 1; i < approxEffResistanceGlobal.size(); ++i)
            approxEffResistanceGlobal[0][u] += approxEffResistanceGlobal[i][u];
        currentRoundResistanceApprox[u] = approxEffResistanceGlobal[0][u];
        resistanceToRoot[u] = (1. - w) * resistanceToRoot[u] + w * currentRoundResistanceApprox[u];
        diagonal[u] = resistanceToRoot[u] - rootCol[root] + 2. * rootCol[u];
    });

    diagonal[root] = rootCol[root];
    diagonal[a] = lpinvColA[a];
    diagonal[b] = lpinvColB[b];
    const double trace = std::accumulate(diagonal.begin(), diagonal.end(), 0.);

    G.parallelForNodes(
        [&](node u) { scoreData[u] = (n_double - 1.) / (n_double * diagonal[u] + trace); });
}

void DynApproxElectricalCloseness::edgeRemoved(node a, node b) {
    assert(!G.hasEdge(a, b));
    assert(hasRun);
    assert(bfsParent[a] != b && bfsParent[b] != a);

    const count n = G.numberOfNodes();
    const double n_double = static_cast<double>(n);

    // precond: this->laplacian is the laplacian of G' = G with edge (a,b)
    //          this->G is the graph without edge (a,b)

    // we need the following values:
    //   r_G'(u,v) for all v, u=pivot - stored in this->resistanceToRoot
    //   r_G'(a,b) - computed via columns/solver. requires old laplacian
    //   E[F | e in t'] - computed via USTs with edge (a,b)
    //   we need the u'th column of the pseudoinverse of L(G) for the diag formula

    // Compute lpinv columns for a and b exactly. (for G', old laplacian)

    const auto [w, lpinvColA, lpinvColB] = computeColumnsExact(a, b);

    // compute u'th column of laplacian pseudoinverse for G

    // update laplacian to G (remove edge (a,b))
    laplacian.setValue(a, a, laplacian(a, a) - 1.);
    laplacian.setValue(b, b, laplacian(b, b) - 1.);
    laplacian.setValue(a, b, laplacian(a, b) + 1.);
    laplacian.setValue(b, a, laplacian(b, a) + 1.);

    solveRootAndResample(a, b, w);

    // Aggregate results
    Vector currentRoundResistanceApprox(n);

    G.parallelForNodes([&](const node u) {
        // Accumulate all results on first thread vector
        for (count i = 1; i < approxEffResistanceGlobal.size(); ++i)
            approxEffResistanceGlobal[0][u] += approxEffResistanceGlobal[i][u];
        currentRoundResistanceApprox[u] = approxEffResistanceGlobal[0][u]; // E[F | e in t']
        resistanceToRoot[u] = // flipped formula compared to edgeAdded
            (resistanceToRoot[u] - w * currentRoundResistanceApprox[u]) / (1. - w);
        diagonal[u] = resistanceToRoot[u] - rootCol[root] + 2. * rootCol[u];
    });

    diagonal[root] = rootCol[root];
    const double trace = std::accumulate(diagonal.begin(), diagonal.end(), 0.);

    G.parallelForNodes(
        [&](node u) { scoreData[u] = (n_double - 1.) / (n_double * diagonal[u] + trace); });
}

std::tuple<double, Vector, Vector> DynApproxElectricalCloseness::computeColumnsExact(node a,
                                                                                     node b) {
    // Compute lpinv columns for a and b exactly.
    const count n = G.numberOfNodes();
    const double n_double = static_cast<double>(n);

    ConjugateGradient<CSRMatrix, DiagonalPreconditioner> cg(tol);
    cg.setupConnected(laplacian);

    Vector rhs_a(n), rhs_b(n), lpinvColA(n), lpinvColB(n);
    G.forNodes([&](node u) {
        rhs_a[u] = -1.0 / n_double;
        rhs_b[u] = -1.0 / n_double;
    });
    rhs_a[a] += 1.;
    rhs_b[b] += 1.;
    auto status = cg.solve(rhs_a, lpinvColA);
    assert(status.converged);
    status = cg.solve(rhs_b, lpinvColB);
    assert(status.converged);

    double sum_a = 0., sum_b = 0.;
    G.forNodes([&](node u) {
        sum_a += lpinvColA[u];
        sum_b += lpinvColB[u];
    });
    lpinvColA -= sum_a / n_double;
    lpinvColB -= sum_b / n_double;

    // update weights
    double w = lpinvColA[a] + lpinvColB[b] - 2. * lpinvColA[b];
    assert(0 <= w && w <= 1.);

    return {w, lpinvColA, lpinvColB};
}

void DynApproxElectricalCloseness::solveRootAndResample(node a, node b, double w) {
    const count n = G.numberOfNodes();
    const double n_double = static_cast<double>(n);

    ConjugateGradient<CSRMatrix, DiagonalPreconditioner> cg(tol);
    cg.setupConnected(laplacian);

    Vector rhs(n);
    G.forNodes([&](node u) { rhs[u] = -1.0 / n_double; });
    rhs[root] += 1.;
    Vector exactRootCol(n);
    cg.solve(rhs, exactRootCol);
    double sum = 0.;
    G.forNodes([&](node u) { sum += exactRootCol[u]; });
    exactRootCol -= sum / n_double;
    rootCol = exactRootCol;

    // Sample and aggregate USTs
    count ustsCurrentRound = std::ceil(w * numberOfUSTs);
    INFO("USTS: ", ustsCurrentRound);

    for (auto &approxEffResistanceLocal : approxEffResistanceGlobal) {
        std::fill(approxEffResistanceLocal.begin(), approxEffResistanceLocal.end(), 0.);
    }

    // update degDist
    degDist[a] = std::uniform_int_distribution<index>(0, G.degree(a) - 1);
    degDist[b] = std::uniform_int_distribution<index>(0, G.degree(b) - 1);

#pragma omp parallel for
    for (omp_index i = 0; i < ustsCurrentRound; ++i) {
        sampleUSTWithEdge(a,
                          b); // for removed edges: edge (a,b) is not in G, but the sample function
                              // does not check this
                              // -> the UST is a UST in G' which is what we need here.
        dfsUST();
        aggregateUST(1. / static_cast<double>(ustsCurrentRound));
    }
}

void DynApproxElectricalCloseness::update(GraphEvent e) {
    if (e.type == GraphEvent::EDGE_ADDITION)
        edgeAdded(e.u, e.v);
    else if (e.type == GraphEvent::EDGE_REMOVAL) {
        // the deleted edge may be part of the bfsTree. If this is the case, we need a new one.
        // Since all our computations depend on this exact tree (via aggregateUST), we need to
        // re-compute some results from previous iterations (either reset the state and run() again,
        // or store all USTs and aggregate them again). for now, we throw an error in this case -
        // current planned use cases do not delete edges from the bfs tree.
        if (bfsParent[e.u] == e.v || bfsParent[e.v] == e.u)
            throw std::runtime_error("Error: Edge removal where the edge removed is in the "
                                     "bfsTree is not supported.");
        else
            edgeRemoved(e.u, e.v);
    } else
        throw std::runtime_error(
            "Error: Events other than EDGE_ADDITION and EDGE_REMOVAL are not supported.");
}

void DynApproxElectricalCloseness::run() {
    // Preprocessing
    computeNodeSequence(pivot);
    computeBFSTree();

    numberOfUSTs = computeNumberOfUSTs();
    rootCol = Vector(G.numberOfNodes());

    solveColumnAndSampleUSTs();
    aggregateResults();

    hasRun = true;
}

} // namespace NetworKit
