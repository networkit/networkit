/*
 * DynApproxElectricalCloseness.cpp
 *
 *  Created on: 17.10.2019
 *     Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *              Alexander van der Grinten <avdgrinten@hu-berlin.de>
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
    : Centrality(G), epsilon(epsilon), delta(delta), kappa(kappa), pivot(pivot),
      bcc(BiconnectedComponents(G)) {

    if (G.isDirected())
        throw std::runtime_error("Error: the input graph must be undirected.");

    if (G.isWeighted())
        throw std::runtime_error("Error: the input graph must be unweighted.");

    if (G.numberOfNodes() < 2)
        throw std::runtime_error("Error: the graph should have at leasts two vertices");

    if (G.upperNodeIdBound() != G.numberOfNodes())
        throw std::runtime_error("Error, graph is not compact. Use getCompactedGraph in GraphTools "
                                 "before using this algorithm.");

    const count n = G.upperNodeIdBound(), threads = omp_get_max_threads();
    statusGlobal.resize(threads, std::vector<NodeStatus>(n, NodeStatus::NOT_VISITED));
    approxEffResistanceGlobal.resize(threads, std::vector<double>(n));
    scoreData.resize(n);
    diagonal.resize(n);
    resistanceToRoot.resize(n);
    generators.reserve(threads);
    for (omp_index i = 0; i < threads; ++i)
        generators.emplace_back(Aux::Random::getURNG());

    degDist.reserve(G.upperNodeIdBound());
    G.forNodes([&](node u) { degDist.emplace_back(0, G.degree(u) - 1); });

    lcaGlobal.resize(threads, std::vector<node>(n, none));
}

DynApproxElectricalCloseness::DynApproxElectricalCloseness(const Graph &G, double epsilon,
                                                           double kappa, node pivot)
    : DynApproxElectricalCloseness(G, epsilon, kappa, pivot,
                                   1.0 / static_cast<double>(G.numberOfNodes())) {}

count DynApproxElectricalCloseness::computeNumberOfUSTs() const {
    return rootEcc * rootEcc
           * static_cast<count>(
               std::ceil(std::log(2.0 * static_cast<double>(G.numberOfEdges()) / delta)
                         / (2.0 * epsilon * epsilon * (1. - kappa) * (1. - kappa))));
}

node DynApproxElectricalCloseness::approxMinEccNode() {
    auto &status = statusGlobal[0];
    std::vector<count> distance(G.upperNodeIdBound());
    std::vector<count> eccLowerBound(G.upperNodeIdBound());

    auto maxDegreeNode = [&]() -> node {
        node maxDegNode = 0;
        count maxDeg = 0;
        G.forNodes([&](const node u) {
            const auto degU = G.degree(u);
            if (degU > maxDeg) {
                maxDeg = degU;
                maxDegNode = u;
            }
        });
        return maxDegNode;
    };

    auto doBFS = [&](const node source) -> node {
        std::queue<node> q;
        q.push(source);
        status[source] = NodeStatus::VISITED;
        distance[source] = 0;
        node farthest = 0;

        do {
            const node u = q.front();
            q.pop();
            eccLowerBound[u] = std::max(eccLowerBound[u], distance[u]);
            farthest = u;

            G.forNeighborsOf(u, [&](const node v) {
                if (status[v] == NodeStatus::NOT_VISITED) {
                    q.push(v);
                    status[v] = NodeStatus::VISITED;
                    distance[v] = distance[u] + 1;
                }
            });

        } while (!q.empty());

        std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);
        return farthest;
    };

    node source = maxDegreeNode();

    for (uint32_t i = 0; i < sweeps; ++i)
        source = doBFS(source);

    // Return node with minimum ecc lower bound
    return static_cast<node>(std::min_element(eccLowerBound.begin(), eccLowerBound.end())
                             - eccLowerBound.begin());
}

void DynApproxElectricalCloseness::computeNodeSequence() {
    // We use thread 0's vector
    auto &status = statusGlobal[0];

    // Compute the biconnected components
    bcc.run();
    auto components = bcc.getComponents();

    std::queue<node> queue;
    std::vector<node> curSequence;

    for (const auto &curComponent : components) {
        // Biconnected components with 2 vertices can be handled trivially.
        if (curComponent.size() == 2) {
            sequences.push_back({curComponent[0], curComponent[1]});
            continue;
        }

        // We take the node with highest degree in the component as source
        node source = curComponent[0];
        for (node u : curComponent) {
            source = G.degree(u) > G.degree(source) ? u : source;
            status[u] = NodeStatus::VISITED;
        }

        // Start a BFS from source to determine the order of the nodes.
        // Later, we need to re-explore the graph again, so in the beginning we mark all nodes to be
        // visited as VISITED, and we set them as NOT_VISITED during the BFS. This avoids a call to
        // std::fill.
        queue.push(source);
        status[source] = NodeStatus::NOT_VISITED;

        do {
            const node u = queue.front();
            queue.pop();
            curSequence.push_back(u);
            G.forNeighborsOf(u, [&](node v) {
                if (status[v] == NodeStatus::VISITED) {
                    status[v] = NodeStatus::NOT_VISITED;
                    queue.push(v);
                }
            });
        } while (!queue.empty());

        sequences.push_back(std::move(curSequence));
        curSequence.clear();
    }

    // Set the root to the highest degree node within the largest biconnected component
    if (pivot == none)
        root = approxMinEccNode();
    else
        root = pivot;

    biAnchor.resize(bcc.numberOfComponents(), none);
    biParent.resize(bcc.numberOfComponents(), none);

#ifdef NETWORKIT_SANITY_CHECKS
    G.forNodes([&](node u) { assert(status[u] == NodeStatus::NOT_VISITED); });
#endif

    // Topological order of biconnected components: tree of biconnected components starting from the
    // root's biconnected component. If the root is in multiple biconnected components, arbitrarily
    // select one of them.
    std::queue<std::pair<node, index>> q;
    const auto &rootComps = bcc.getComponentsOfNode(root);
    q.push({root, *(rootComps.begin())});

    topOrder.reserve(bcc.numberOfComponents());
    topOrder.insert(topOrder.begin(), rootComps.begin(), rootComps.end());

    std::vector<count> distance(G.upperNodeIdBound());
    status[root] = NodeStatus::VISITED;
    rootEcc = 0;

    do {
        const auto front = q.front();
        q.pop();
        G.forNeighborsOf(front.first, [&](const node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                distance[v] = distance[front.first] + 1;
                rootEcc = std::max(rootEcc, distance[v]);
                const auto &vComps = bcc.getComponentsOfNode(v);
                for (const node vComponentIndex : vComps) {
                    // Check if a new biconnected components has been found.
                    if (vComponentIndex != front.second && biAnchor[vComponentIndex] == none
                        && rootComps.find(vComponentIndex) == rootComps.end()) {
                        // The anchor cannot be the root, because the anchor does not have a parent.
                        // We handle biAnchor = none cases later.
                        biAnchor[vComponentIndex] = (v == root) ? none : v;
                        biParent[vComponentIndex] = front.second;
                        topOrder.push_back(vComponentIndex);
                    }
                }
                q.push({v, *(vComps.begin())});
                status[v] = NodeStatus::VISITED;
            }
        });
    } while (!q.empty());

#ifdef NETWORKIT_SANITY_CHECKS
    G.forNodes([&](node u) { assert(status[u] == NodeStatus::VISITED); });
    assert(topOrder.size() == bcc.numberOfComponents());
    assert(std::unordered_set<index>(topOrder.begin(), topOrder.end()).size()
           == bcc.numberOfComponents());
#endif
}

void DynApproxElectricalCloseness::computeBFSTree() {
    // Using thread 0's vector
    auto &status = statusGlobal[0];
    auto n = G.numberOfNodes();
    bfsTree.parent.resize(n);
    bfsTree.child.resize(n);
    bfsTree.sibling.resize(n);
    std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);
    std::fill(bfsTree.parent.begin(), bfsTree.parent.end(), none);
    std::fill(bfsTree.child.begin(), bfsTree.child.end(), none);
    std::fill(bfsTree.sibling.begin(), bfsTree.sibling.end(), none);

    std::queue<node> queue;
    queue.push(root);
    status[root] = NodeStatus::VISITED;

    do {
        const node currentNode = queue.front();
        queue.pop();
        node previous = none;
        node first = none;
        G.forNeighborsOf(currentNode, [&](const node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                status[v] = NodeStatus::VISITED;
                queue.push(v);
                bfsTree.parent[v] = currentNode;
                if (previous != none) {
                    bfsTree.sibling[previous] = v;
                }
                previous = v;
                if (first == none) {
                    bfsTree.child[currentNode] = v;
                    first = v;
                }
            }
        });
    } while (!queue.empty());

#ifdef NETWORKIT_SANITY_CHECKS
    checkBFSTree();
#endif
}

void DynApproxElectricalCloseness::sampleUST(Tree &result) {
    // Getting thread-local vectors
    auto &status = statusGlobal[omp_get_thread_num()];
    auto &parent = result.parent;
    auto &childPtr = result.child;
    auto &siblingPtr = result.sibling;

    auto n = G.numberOfNodes();
    parent.resize(n);
    childPtr.resize(n);
    siblingPtr.resize(n);
    std::fill(status.begin(), status.end(), NodeStatus::NOT_IN_COMPONENT);
    std::fill(parent.begin(), parent.end(), none);
    std::fill(childPtr.begin(), childPtr.end(), none);
    std::fill(siblingPtr.begin(), siblingPtr.end(), none);

    auto &generator = generators[omp_get_thread_num()];

    // Iterate over the biconnected components in their topological order.
    for (const auto componentIndex : topOrder) {
        // Current component, sorted by vertex degree.
        const auto &sequence = sequences[componentIndex];
        auto curAnchor = biAnchor[componentIndex];

        // Finds the parent of the current anchor node i.e., the anchor's neighbor that is in the
        // parent component.
        auto updateParentOfAnchor = [&]() -> void {
            for (const node v : G.neighborRange(curAnchor)) {
                const auto &vComps = bcc.getComponentsOfNode(v);
                if (vComps.find(biParent[componentIndex]) != vComps.end()) {
                    parent[curAnchor] = v;
                    break;
                }
            }
            assert(parent[curAnchor] != none);
        };

        if (sequence.size() == 2) {
            // Happens when the current component is the root component in the topological
            // order. In this case, the root plays the anchor's role.
            if (curAnchor == none) {
                const node v = (sequence[0] == root) ? sequence[1] : sequence[0];
                assert(sequence.front() == root || sequence.back() == root);
                assert(v != root);
                parent[v] = root;
            } else {
                const node v = (sequence[0] == curAnchor) ? sequence[1] : sequence[0];
                assert(v != curAnchor);
                assert(v != root);
                parent[v] = curAnchor;
                if (parent[curAnchor] == none)
                    updateParentOfAnchor();
            }
#ifdef NETWORKIT_SANITY_CHECKS
            checkTwoNodesSequence(sequence, parent);
#endif
            continue;
        }

        // We start building the spanning tree from the first node of the
        // sequence.
        // Root of the spanning tree
        status[sequence[0]] = NodeStatus::IN_SPANNING_TREE;
        const node curAnchorParent = (curAnchor != none) ? parent[curAnchor] : none;

        // All the remaining nodes in the components need to be visited.
        std::for_each(sequence.begin() + 1, sequence.end(),
                      [&status](node u) { status[u] = NodeStatus::NOT_VISITED; });

        count nodesInSpanningTree = 1;
        // Iterate over the remaining nodes to create the spanning tree.
        for (auto it = sequence.begin() + 1; it != sequence.end(); ++it) {
            const node walkStart = *it;
            if (status[walkStart] == NodeStatus::IN_SPANNING_TREE)
                // Node already added to the spanning tree
                continue;

            node currentNode = walkStart;

            // Start a new random walk from the current node.
            do {
                // Get a random neighbor within the component
                node randomNeighbor;
                do {
                    randomNeighbor = G.getIthNeighbor(currentNode, degDist[currentNode](generator));
                } while (status[randomNeighbor] == NodeStatus::NOT_IN_COMPONENT);

                assert(randomNeighbor != none);
                parent[currentNode] = randomNeighbor;
                currentNode = randomNeighbor;

            } while (status[currentNode] != NodeStatus::IN_SPANNING_TREE);

            // Last node encountered in the random walk (it is in the spanning tree);
            const node walkEnd = currentNode;
            assert(status[walkEnd] == NodeStatus::IN_SPANNING_TREE);
            // Add the random walk to the spanning tree; eventually, reverse the path if the
            // anchor/root is encountered
            for (currentNode = walkStart; currentNode != walkEnd;
                 currentNode = parent[currentNode]) {

                status[currentNode] = NodeStatus::IN_SPANNING_TREE;
                ++nodesInSpanningTree;
                if (currentNode == curAnchor || currentNode == root) {

                    // Anchor of current component in the walk, we have to reverse the
                    // parent pointers
                    node next = parent[currentNode];
                    node nextParent = none;
                    node tmp = currentNode;
                    do {
                        status[next] = NodeStatus::IN_SPANNING_TREE;
                        nextParent = parent[next];
                        parent[next] = currentNode;
                        currentNode = next;
                        next = nextParent;
                    } while (next != none);

                    if (tmp == root)
                        parent[root] = none;
                    else if (parent[curAnchor] == none)
                        // Not none if articulation node visited by the parent component before
                        updateParentOfAnchor();
                    else
                        parent[curAnchor] = curAnchorParent;

                    break;
                }
            }

            if (nodesInSpanningTree == sequence.size())
                break;
        }

        for (node u : sequence)
            status[u] = NodeStatus::NOT_IN_COMPONENT;
    }

    count visitedNodes = 0;
    for (node u : G.nodeRange()) {
        while (status[u] == NodeStatus::NOT_IN_COMPONENT) {
            status[u] = NodeStatus::NOT_VISITED;
            ++visitedNodes;
            node parentU = parent[u];
            if (parent[u] != none) {
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
    checkUST(result);
#endif
}

void DynApproxElectricalCloseness::sampleUSTWithEdge(Tree &result, node a, node b) {
    auto n = G.numberOfNodes();
    assert(a < n && b < n);

    auto &parent = result.parent;
    auto &status = statusGlobal[omp_get_thread_num()];
    auto &childPtr = result.child;
    auto &siblingPtr = result.sibling;

    parent.resize(n, none);
    childPtr.resize(n, none);
    siblingPtr.resize(n, none);
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
    checkUST(result);
#endif
}

void DynApproxElectricalCloseness::setDFSTimes(Tree &tree) {
    if (tree.timesComputed) {
        return;
    }
    tree.tVisit.resize(G.numberOfNodes());
    tree.tFinish.resize(G.numberOfNodes());

    std::stack<std::pair<node, node>> stack;
    stack.push({root, tree.child[root]});
    count timestamp = 0;

    do {
        // v is a child of u that has not been visited yet.
        const node u = stack.top().first;
        const node v = stack.top().second;

        if (v == none) {
            stack.pop();
            tree.tFinish[u] = ++timestamp;
        } else {
            stack.top().second = tree.sibling[v];
            tree.tVisit[v] = ++timestamp;
            stack.push({v, tree.child[v]});
        }
    } while (!stack.empty());

    tree.timesComputed = true;
}

void DynApproxElectricalCloseness::aggregateUST(Tree &tree, double weight) {
#ifdef NETWORKIT_SANITY_CHECKS
    checkTimeStamps(tree);
#endif

    auto &approxEffResistance = approxEffResistanceGlobal[omp_get_thread_num()];
    const auto &tVisit = tree.tVisit;
    const auto &tFinish = tree.tFinish;
    const auto &parent = tree.parent;

    // Doing aggregation
    G.forNodes([&](const node u) {
        node p = bfsTree.parent[u];
        node c = u;

        auto goUp = [&]() -> void {
            c = p;
            p = bfsTree.parent[p];
        };

        while (p != none) {
            // Edge in BSTree: e1 -> e2
            node e1 = p, e2 = c;
            bool reverse = false;
            if (e1 != parent[e2]) {
                if (e2 != parent[e1]) {
                    goUp();
                    continue;
                }
                std::swap(e1, e2);
                reverse = true;
            }

            if (tVisit[u] >= tVisit[e2] && tFinish[u] <= tFinish[e2])
                approxEffResistance[u] += reverse ? -weight : weight;

            goUp();
        }
    });
}

void DynApproxElectricalCloseness::edgeAdded(node a, node b) {
    assert(G.hasEdge(a, b));
    assert(hasRun);

    // Compute lpinv columns for a and b exactly.
    const count n = G.numberOfNodes();
    const double n_double = static_cast<double>(n);

    /// double laa = laplacian(a, a), lab = laplacian(a, b), lba = laplacian(b, a), lbb =
    /// laplacian(b, b);
    laplacian.setValue(a, a, laplacian(a, a) + 1.);
    laplacian.setValue(b, b, laplacian(b, b) + 1.);
    laplacian.setValue(a, b, laplacian(a, b) - 1.);
    laplacian.setValue(b, a, laplacian(b, a) - 1.);

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
        Tree tree;
        sampleUSTWithEdge(tree, a, b);
        setDFSTimes(tree);
        aggregateUST(tree, 1. / static_cast<double>(ustsCurrentRound));
    }

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
    assert(bfsTree.parent[a] != b && bfsTree.parent[b] != a);

    // precond: this->laplacian is the laplacian of G' = G with edge (a,b)
    //          this->G is the graph without edge (a,b)

    // we need the following values:
    //   r_G'(u,v) for all v, u=pivot - stored in this->resistanceToRoot
    //   r_G'(a,b) - computed via columns/solver. requires old laplacian
    //   E[F | e in t'] - computed via USTs with edge (a,b)
    //   we need the u'th column of the pseudoinverse of L(G) for the diag formula

    // Compute lpinv columns for a and b exactly. (for G', old laplacian)
    const count n = G.numberOfNodes();
    const double n_double = static_cast<double>(n);

    ConjugateGradient<CSRMatrix, DiagonalPreconditioner> cg_old(tol);
    cg_old.setupConnected(laplacian);

    Vector rhs_a(n), rhs_b(n), lpinvColA(n), lpinvColB(n);
    G.forNodes([&](node u) {
        rhs_a[u] = -1.0 / n_double;
        rhs_b[u] = -1.0 / n_double;
    });
    rhs_a[a] += 1.;
    rhs_b[b] += 1.;
    auto status = cg_old.solve(rhs_a, lpinvColA);
    assert(status.converged);
    status = cg_old.solve(rhs_b, lpinvColB);
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

    // compute u'th column of laplacian pseudoinverse for G

    // update laplacian to G (remove edge (a,b))
    laplacian.setValue(a, a, laplacian(a, a) - 1.);
    laplacian.setValue(b, b, laplacian(b, b) - 1.);
    laplacian.setValue(a, b, laplacian(a, b) + 1.);
    laplacian.setValue(b, a, laplacian(b, a) + 1.);

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
        Tree tree;
        sampleUSTWithEdge(tree, a, b);
        setDFSTimes(tree);
        aggregateUST(tree, 1. / static_cast<double>(ustsCurrentRound));
    }

    // Aggregate results
    Vector currentRoundResistanceApprox(n);

    G.parallelForNodes([&](const node u) {
        // Accumulate all results on first thread vector
        for (count i = 1; i < approxEffResistanceGlobal.size(); ++i)
            approxEffResistanceGlobal[0][u] += approxEffResistanceGlobal[i][u];
        currentRoundResistanceApprox[u] = approxEffResistanceGlobal[0][u]; // E[F | e in t']
        resistanceToRoot[u] =
            (resistanceToRoot[u] - w * currentRoundResistanceApprox[u]) / (1. - w);
        diagonal[u] = resistanceToRoot[u] - rootCol[root] + 2. * rootCol[u];
    });

    diagonal[root] = rootCol[root];
    const double trace = std::accumulate(diagonal.begin(), diagonal.end(), 0.);

    G.parallelForNodes(
        [&](node u) { scoreData[u] = (n_double - 1.) / (n_double * diagonal[u] + trace); });
}

void DynApproxElectricalCloseness::update(GraphEvent e) {
    if (e.type == GraphEvent::EDGE_ADDITION)
        edgeAdded(e.u, e.v);
    else if (e.type == GraphEvent::EDGE_REMOVAL) {
        // the deleted edge may be part of the bfsTree. If this is the case, we need a new one.
        // Since all our computations depend on this exact tree (via aggregateUST), we need to
        // re-compute some results from previous iterations.
        // for now, we throw an error in this case - current planned use cases do not delete
        // edges from the bfs tree.
        if (bfsTree.parent[e.u] == e.v || bfsTree.parent[e.v] == e.u)
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
    computeNodeSequence();
    computeBFSTree();

    numberOfUSTs = computeNumberOfUSTs();
    rootCol = Vector(G.numberOfNodes());

    double weight = 1. / static_cast<double>(numberOfUSTs);

#pragma omp parallel
    {
        // Thread 0 solves the linear system
        if (omp_get_thread_num() == 0) {
            const count n = G.numberOfNodes();

            laplacian = CSRMatrix::laplacianMatrix(G);
            Diameter diamAlgo(G, DiameterAlgo::ESTIMATED_RANGE, 0);
            diamAlgo.run();
            // Getting diameter upper bound
            const auto diam = static_cast<double>(diamAlgo.getDiameter().second);
            tol = epsilon * kappa
                  / (std::sqrt(static_cast<double>((n * G.numberOfEdges())) * std::log(n)) * diam
                     * 3.);
            ConjugateGradient<CSRMatrix, DiagonalPreconditioner> cg(tol);
            cg.setupConnected(laplacian);

            Vector rhs(n);
            G.forNodes([&](node u) { rhs[u] = -1.0 / static_cast<double>(n); });
            rhs[root] += 1.;
            cg.solve(rhs, rootCol);

            double sum = 0.0;
            G.forNodes([&](node u) { sum += rootCol[u]; });
            rootCol -= sum / static_cast<double>(n);
        }

        // All threads sample and aggregate USTs in parallel
#pragma omp for
        for (omp_index i = 0; i < numberOfUSTs; ++i) {
            Tree tree;
            sampleUST(tree);
            setDFSTimes(tree);
            aggregateUST(tree, weight);
        }
    }

    // Aggregating thread-local results
    G.parallelForNodes([&](const node u) {
        // Accumulate all results on thread 0 vector
        for (count i = 1; i < approxEffResistanceGlobal.size(); ++i)
            approxEffResistanceGlobal[0][u] += approxEffResistanceGlobal[i][u];
        resistanceToRoot[u] = approxEffResistanceGlobal[0][u];
        diagonal[u] = resistanceToRoot[u] - rootCol[root] + 2. * rootCol[u];
    });

    diagonal[root] = rootCol[root];
    const double trace = std::accumulate(diagonal.begin(), diagonal.end(), 0.);
    const double n = G.numberOfNodes();

    G.parallelForNodes([&](node u) { scoreData[u] = (n - 1.) / (n * diagonal[u] + trace); });

    hasRun = true;
}

#ifdef NETWORKIT_SANITY_CHECKS
/*
 * Methods for sanity check.
 */
void DynApproxElectricalCloseness::checkUST(const Tree &tree) const {
    std::vector<bool> visitedNodes(G.upperNodeIdBound());
    const auto &parent = tree.parent;
    const auto &childPtr = tree.child;
    const auto &siblingPtr = tree.sibling;

    // To debug
    G.forNodes([&](node u) {
        if (childPtr[u] != none) {
            assert(u == parent[childPtr[u]]);
        }
        if (siblingPtr[u] != none) {
            assert(parent[u] == parent[siblingPtr[u]]);
        }
        if (u == root) {
            assert(parent[u] == none);
        } else {
            assert(parent[u] != none);
            std::fill(visitedNodes.begin(), visitedNodes.end(), 0);
            visitedNodes[u] = 1;
            do {
                u = parent[u];
                assert(!visitedNodes[u]);
                visitedNodes[u] = 1;
            } while (u != root);
        }
    });
}

void DynApproxElectricalCloseness::checkBFSTree() const {
    G.forNodes([&](node u) {
        if (u == root) {
            assert(bfsTree.parent[u] == none);
        } else {
            std::vector<bool> visited(G.upperNodeIdBound());
            do {
                assert(!visited[u]);
                visited[u] = true;
                assert(bfsTree.parent[u] != none);
                u = bfsTree.parent[u];
            } while (u != root);
        }
    });
}

void DynApproxElectricalCloseness::checkTwoNodesSequence(const std::vector<node> &sequence,
                                                         std::vector<node> &parent) const {
    for (node u : sequence) {
        if (u == root) {
            assert(parent[u] == none);
        } else {
            std::vector<bool> visited(G.upperNodeIdBound());
            visited[u] = true;
            do {
                u = parent[u];
                assert(!visited[u]);
                visited[u] = true;
            } while (u != root);
        }
    }
}

void DynApproxElectricalCloseness::checkTimeStamps(const Tree &tree) const {
    const auto &tVisit = tree.tVisit;
    const auto &tFinish = tree.tFinish;
    G.forNodes([&](const node u) {
        assert(tVisit[u] < tFinish[u]);
        if (u == root)
            assert(tVisit[u] == 0);
        else
            assert(tVisit[u] > 0);
        assert(tVisit[u] < 2 * G.numberOfNodes());
    });
}

#endif // NETWORKIT_SANITY_CHECKS

} // namespace NetworKit
