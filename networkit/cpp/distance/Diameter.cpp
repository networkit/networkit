/*
 * Diameter.cpp
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <span>
#include <utility>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/structures/Partition.hpp>

#include "VertexDiameterBound.hpp"

namespace NetworKit {

namespace {

/*
 * Bounds on the eccentricities of the nodes of a directed graph that does not have to be strongly
 * connected. The forward eccentricity of a node is its largest finite distance to another node, and
 * the backward eccentricity the largest finite distance from another node. The diameter is the
 * largest eccentricity of either kind, so it is enough to bound the eccentricities of one kind.
 *
 * This extends the idea of "On computing the diameter of real-world directed (weighted) graphs"
 * by Crescenzi et al. (SEA 2012), which assumes a strongly connected graph, in two ways: bounds
 * follow the edges between strongly connected components, and a longest shortest path that exceeds
 * the lower bound has to start and end at nodes whose eccentricity can exceed it.
 */
class DirectedEccentricityBounds {
public:
    enum Direction { FORWARD = 0, BACKWARD = 1 };

    static constexpr Direction opposite(Direction direction) {
        return direction == FORWARD ? BACKWARD : FORWARD;
    }

    explicit DirectedEccentricityBounds(const Graph &G)
        : G(G), transposed(GraphTools::transpose(G)), scc(G) {
        scc.run();
        membersOfComponent.resize(scc.numberOfComponents());
        G.forNodes(
            [&](node u) -> void { membersOfComponent[scc.componentOfNode(u)].push_back(u); });
        bounds[FORWARD].assign(G.upperNodeIdBound(), G.numberOfNodes() - 1);
        bounds[BACKWARD] = bounds[FORWARD];

        const auto largestComponent = std::ranges::max_element(
            membersOfComponent, {},
            [](std::span<const node> members) -> count { return members.size(); });
        pivotComponent = static_cast<index>(largestComponent - membersOfComponent.begin());
        const std::span<const node> members = membersOfComponent[pivotComponent];
        const node pivot =
            *std::ranges::max_element(members, {}, [&](node u) -> count { return totalDegree(u); });
        distancesFromPivot[FORWARD] = sweep(FORWARD, pivot);
        distancesFromPivot[BACKWARD] = sweep(BACKWARD, pivot);
        propagate();
    }

    count lowerBound() const { return lower; }

    // Runs a BFS from source in the given direction and returns the distances, none for the
    // unreachable nodes.
    std::vector<count> sweep(Direction direction, node source) {
        std::vector<count> distances(G.upperNodeIdBound(), none);
        count eccentricity = 0;
        Traversal::BFSfrom(graph(direction), source, [&](node u, count dist) -> void {
            distances[u] = dist;
            eccentricity = std::max(eccentricity, dist);
        });
        lower = std::max(lower, eccentricity);
        bounds[direction][source] = std::min(bounds[direction][source], eccentricity);

        // All the nodes of a component are connected to the same nodes, so they can go via source.
        const Direction otherDirection = opposite(direction);
        const count eccentricityOfSource = bounds[otherDirection][source];
        for (node u : membersOfComponent[scc.componentOfNode(source)]) {
            bounds[otherDirection][u] =
                std::min(bounds[otherDirection][u], distances[u] + eccentricityOfSource);
        }
        return distances;
    }

    // Uses that a pair of nodes that is further apart than target consists of candidates.
    void refine(count target) {
        for (Direction direction : {FORWARD, BACKWARD}) {
            const Direction otherDirection = opposite(direction);
            count furthestCandidate = 0;
            G.forNodes([&](node u) -> void {
                if (bounds[otherDirection][u] > target
                    && distancesFromPivot[direction][u] != none) {
                    furthestCandidate =
                        std::max(furthestCandidate, distancesFromPivot[direction][u]);
                }
            });
            for (node u : membersOfComponent[pivotComponent]) {
                bounds[direction][u] = std::min(
                    bounds[direction][u],
                    std::max(target, distancesFromPivot[otherDirection][u] + furthestCandidate));
            }
        }
        propagate();
    }

    count numberOfCandidates(Direction direction, count target) const {
        count candidates = 0;
        G.forNodes([&](node u) -> void { candidates += bounds[direction][u] > target; });
        return candidates;
    }

    count largestBound(Direction direction) const {
        count largest = 0;
        G.forNodes([&](node u) -> void { largest = std::max(largest, bounds[direction][u]); });
        return largest;
    }

    node nodeWithLargestBound(Direction direction) const {
        node best = none;
        G.forNodes([&](node u) -> void {
            if (best == none
                || std::make_pair(bounds[direction][u], totalDegree(u))
                       > std::make_pair(bounds[direction][best], totalDegree(best))) {
                best = u;
            }
        });
        return best;
    }

private:
    const Graph &graph(Direction direction) const { return direction == FORWARD ? G : transposed; }

    count totalDegree(node u) const { return G.degreeOut(u) + G.degreeIn(u); }

    // A shortest path visits each component at most once, so a node is bounded by the path through
    // the rest of its component plus the bound of the node reached by the edge that leaves it.
    void propagateComponent(Direction direction, std::span<const node> members) {
        count boundAfterLeaving = 0;
        for (node u : members) {
            graph(direction).forNeighborsOf(u, [&](node v) -> void {
                if (scc.componentOfNode(u) != scc.componentOfNode(v)) {
                    boundAfterLeaving = std::max(boundAfterLeaving, bounds[direction][v] + 1);
                }
            });
        }
        const count bound = members.size() - 1 + boundAfterLeaving;
        for (node u : members) {
            bounds[direction][u] = std::min(bounds[direction][u], bound);
        }
    }

    // Components are numbered in reverse topological order, so the components that a component has
    // edges to are processed before it.
    void propagate() {
        for (auto it = membersOfComponent.begin(); it != membersOfComponent.end(); ++it) {
            propagateComponent(FORWARD, *it);
        }
        for (auto it = membersOfComponent.rbegin(); it != membersOfComponent.rend(); ++it) {
            propagateComponent(BACKWARD, *it);
        }
    }

    const Graph &G;
    const Graph transposed;
    StronglyConnectedComponents scc;
    std::vector<std::vector<node>> membersOfComponent;
    std::array<std::vector<count>, 2> bounds;
    std::array<std::vector<count>, 2> distancesFromPivot;
    index pivotComponent = 0;
    count lower = 0;
};

} // namespace

Diameter::Diameter(const Graph &G, DiameterAlgo algo, double error, count nSamples)
    : Algorithm(), G(&G), error(error), nSamples(nSamples) {
    if (algo == DiameterAlgo::AUTOMATIC) {
        this->algo = DiameterAlgo::exact;
    } else {
        this->algo = algo;
        if (this->algo == DiameterAlgo::ESTIMATED_RANGE) {
            if (error == -1.f)
                throw std::invalid_argument(
                    "For Diameter::estimatedRange the parameter error(>=0) has to be supplied");
        } else if (this->algo == DiameterAlgo::ESTIMATED_SAMPLES) {
            if (nSamples == 0)
                throw std::invalid_argument(
                    "For Diameter::estimatedSamples the parameter nSamples(>0) has to be supplied");
        }
    }
}

void Diameter::run() {
    diameterBounds = {0, 0};

    if (G->isEmpty()) {
        hasRun = true;
        return;
    }

    if (algo == DiameterAlgo::EXACT) {
        std::get<0>(diameterBounds) = this->exactDiameter(*G);
    } else if (algo == DiameterAlgo::ESTIMATED_RANGE) {
        diameterBounds = this->estimatedDiameterRange(*G, error);
    } else if (algo == DiameterAlgo::ESTIMATED_SAMPLES) {
        std::get<0>(diameterBounds) = this->estimatedVertexDiameter(*G, nSamples);
    } else if (algo == DiameterAlgo::ESTIMATED_PEDANTIC) {
        std::get<0>(diameterBounds) = this->estimatedVertexDiameterPedantic(*G);
    } else {
        throw std::runtime_error("should never reach this code as the algorithm should be set "
                                 "correctly in the constructor or fail there");
    }
    hasRun = true;
}

std::pair<count, count> Diameter::getDiameter() const {
    return diameterBounds;
}

edgeweight Diameter::exactDiameter(const Graph &G) {
    using namespace std;

    Aux::SignalHandler handler;

    edgeweight diameter = 0.0;

    if (!G.isWeighted()) {
        std::tie(diameter, std::ignore) = estimatedDiameterRange(G, 0);
    } else {
        G.forNodes([&](node v) {
            handler.assureRunning();
            Dijkstra dijkstra(G, v);
            dijkstra.run();
            auto distances = dijkstra.getDistances();
            G.forNodes([&](node u) {
                if (diameter < distances[u]
                    && (!G.isDirected()
                        || distances[u] != std::numeric_limits<edgeweight>::max())) {
                    diameter = distances[u];
                }
            });
        });
    }

    if (diameter == std::numeric_limits<edgeweight>::max()) {
        throw std::runtime_error("Graph not connected - diameter is infinite");
    }
    return diameter;
}

std::pair<edgeweight, edgeweight> Diameter::estimatedDiameterRange(const Graph &G, double error) {
    if (G.isWeighted()) {
        WARN("The input graph is weighted, but this algorithm ignores weights.");
    }
    if (G.isDirected()) {
        return estimatedDirectedDiameterRange(G, error);
    }

    Aux::SignalHandler handler;
    /*
     * This is an implementation of a slightly modified version of the exactSumSweep algorithm as
     * presented in Fast diameter and radius BFS-based computation in (weakly connected) real-world
     * graphs: With an application to the six degrees of separation games by Michele Borassi,
     * Pierluigi Crescenzi, Michel Habib, Walter A. Kosters, Andrea Marino, Frank W. Takes
     * http://www.sciencedirect.com/science/article/pii/S0304397515001644
     */

    std::vector<count> eccLowerBound(G.upperNodeIdBound()), eccUpperBound(G.upperNodeIdBound());
    std::vector<bool> finished(G.upperNodeIdBound());

    G.parallelForNodes([&](node u) { eccUpperBound[u] = G.numberOfNodes(); });

    ConnectedComponents comp(G);
    comp.run();
    count numberOfComponents = comp.numberOfComponents();

    std::vector<count> distFirst;
    std::vector<count> ecc(numberOfComponents, 0);
    std::vector<count> distances(G.upperNodeIdBound(), 0);

    count numBFS = 0;

    auto runBFS = [&](const std::vector<node> &startNodes) {
        ++numBFS;
        distFirst.clear();
        distFirst.resize(numberOfComponents, 0);
        std::vector<bool> foundFirstDeg2Node(numberOfComponents, false);

        Traversal::BFSfrom(G, startNodes.begin(), startNodes.end(), [&](node v, count dist) {
            distances[v] = dist;

            index c = comp.componentOfNode(v);
            ecc[c] = std::max(dist, ecc[c]);

            if (!foundFirstDeg2Node[c] && G.degree(v) > 1) {
                foundFirstDeg2Node[c] = true;
                distFirst[c] = dist;
            }
        });

        G.forNodes([&](node u) {
            if (finished[u])
                return;

            auto c = comp.componentOfNode(u);

            auto eccValue = std::max(distances[u], ecc[c] - distances[u]);
            eccLowerBound[u] = std::max(eccValue, eccLowerBound[u]);

            if (distances[u] <= distFirst[c]) {
                eccUpperBound[u] = eccValue;
            } else {
                eccUpperBound[u] =
                    std::min(eccUpperBound[u], distances[u] + ecc[c] - 2 * distFirst[c]);
            }

            finished[u] = (eccUpperBound[u] == eccLowerBound[u]);
        });

        ecc.clear();
        ecc.resize(numberOfComponents, 0);
    };

    auto diameterBounds = [&]() {
        count maxExact = *std::max_element(eccLowerBound.begin(), eccLowerBound.end());
        count maxPotential = *std::max_element(eccUpperBound.begin(), eccUpperBound.end());
        return std::make_pair(maxExact, maxPotential);
    };

    std::vector<node> startNodes(numberOfComponents, 0), maxDist(numberOfComponents, 0);
    std::vector<count> maxDeg(numberOfComponents, 0);

    count lb = 0, ub = G.numberOfNodes();

    auto printStartNodes = [&]() {
        DEBUG("Start nodes (lb = ", lb, ", ub = ", ub, "): ");
        for (node u : startNodes) {
            (void)u; // prevent unused variable warning
            DEBUG("Node ", u, " with distance ", distances[u], ", lower bound ", eccLowerBound[u],
                  ", upper bound ", eccUpperBound[u]);
        }
    };

    // for each component, find the node with the maximum degree and add it as start node
    G.forNodes([&](node v) {
        count d = G.degree(v);
        count c = comp.componentOfNode(v);
        if (d >= maxDeg[c]) {
            startNodes[c] = v;
            maxDeg[c] = d;
        }
    });

    handler.assureRunning();

    runBFS(startNodes);

    std::tie(lb, ub) = diameterBounds();

    for (index i = 0; i < 2 * G.numberOfNodes() && ub > (lb + error * lb); ++i) {
        handler.assureRunning();
        startNodes.clear();
        startNodes.resize(numberOfComponents, none);

        if ((i % 2) == 0) {
            G.forNodes([&](node u) {
                auto c = comp.componentOfNode(u);
                if (startNodes[c] == none
                    || std::tie(eccUpperBound[u], distances[u])
                           > std::tie(eccUpperBound[startNodes[c]], distances[startNodes[c]])) {
                    startNodes[c] = u;
                }
            });
        } else {
            G.forNodes([&](node u) {
                auto c = comp.componentOfNode(u);
                // Idea: we select a node that is central (i.e. has a low lower bound) but that is
                // also close to the previous, non-central node. More generally, the best upper
                // bound we can hope for a node v is eccLowerBound[u] + distance(u, v). We select
                // the node the provides the best upper bound for the previous node u in the hope
                // that in its neighborhood there are more nodes for which the bounds can be
                // decreased. Among all these nodes we select the one that has the largest distance
                // to the previous start node.
                if (startNodes[c] == none) {
                    startNodes[c] = u;
                } else {
                    auto compU = eccLowerBound[u] + distances[u],
                         compStart = eccLowerBound[startNodes[c]] + distances[startNodes[c]];
                    if (distances[u] > distFirst[c]
                        && (compU < compStart
                            || (compU == compStart && distances[u] > distances[startNodes[c]]))) {
                        startNodes[c] = u;
                    }
                }
            });
        }

        handler.assureRunning();

        printStartNodes();
        runBFS(startNodes);

        std::tie(lb, ub) = diameterBounds();
    }

    INFO(numBFS, " BFS used");

    return {lb, ub};
}

std::pair<edgeweight, edgeweight> Diameter::estimatedDirectedDiameterRange(const Graph &G,
                                                                           double error) {
    using Bounds = DirectedEccentricityBounds;

    Aux::SignalHandler handler;
    Bounds bounds(G);
    // Terminates after at most 2n sweeps: bounds only decrease and the target only increases, so
    // nodes never become candidates again, and each sweep drops its source's bound to its exact
    // eccentricity, which is at most the lower bound, so it stops being a candidate.
    while (true) {
        handler.assureRunning();
        const count lower = bounds.lowerBound();
        const count target =
            std::max(lower, static_cast<count>(static_cast<double>(lower) * (1. + error)));
        bounds.refine(target);

        const count forward = bounds.numberOfCandidates(Bounds::FORWARD, target);
        const count backward = bounds.numberOfCandidates(Bounds::BACKWARD, target);
        const Bounds::Direction direction =
            forward <= backward ? Bounds::FORWARD : Bounds::BACKWARD;
        if (std::min(forward, backward) == 0) {
            return {bounds.lowerBound(), bounds.largestBound(direction)};
        }
        bounds.sweep(direction, bounds.nodeWithLargestBound(direction));
    }
}

edgeweight Diameter::estimatedVertexDiameter(const Graph &G, count samples) {

    // TODO: consider weights

    const Graph transposed = G.isDirected() ? GraphTools::transpose(G) : Graph();

    auto estimateFrom = [&](node v) -> count {
        BFS bfs(G, v);
        bfs.run();
        if (G.isDirected()) {
            // A path through v consists of a path to v and a path from v.
            BFS transposedBfs(transposed, v);
            transposedBfs.run();
            const auto largestDistance = [](std::span<const edgeweight> distances) -> edgeweight {
                edgeweight largest = 0;
                for (edgeweight dist : distances) {
                    if (dist != std::numeric_limits<edgeweight>::max()) {
                        largest = std::max(largest, dist);
                    }
                }
                return largest;
            };
            return static_cast<count>(largestDistance(bfs.getDistances())
                                      + largestDistance(transposedBfs.getDistances()))
                   + 1;
        }

        return DiameterDetails::vertexDiameterBound<edgeweight>(bfs.getDistances());
    };

    edgeweight vdMax = 0;
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(samples); ++i) {
        node u = GraphTools::randomNode(G);
        edgeweight vd = estimateFrom(u);
        DEBUG("sampled vertex diameter from node ", u, ": ", vd);
#pragma omp critical
        vdMax = std::max(vdMax, vd);
    }

    return vdMax;
}

edgeweight Diameter::estimatedVertexDiameterPedantic(const Graph &G) {
    if (G.isDirected()) {
        using Bounds = DirectedEccentricityBounds;
        Bounds bounds(G);
        bounds.refine(bounds.lowerBound());
        return std::min(bounds.largestBound(Bounds::FORWARD), bounds.largestBound(Bounds::BACKWARD))
               + 1;
    }

    count vd = 0;
    if (!G.isWeighted()) {
        std::vector<bool> visited(G.upperNodeIdBound(), false);
        // perform breadth-first searches
        G.forNodes([&](node u) {
            if (visited[u] == false) {
                count maxDist = 0, maxDist2 = 0;
                Traversal::BFSfrom(G, u, [&](node v, count dist) {
                    visited[v] = true;
                    if (dist > maxDist) {
                        maxDist2 = maxDist;
                        maxDist = dist;
                    } else if (dist > maxDist2) {
                        maxDist2 = dist;
                    }
                });
                if (maxDist + maxDist2 > vd) {
                    vd = maxDist + maxDist2;
                }
                assert(visited[u] == true);
            }
        });
        vd++; // we need the number of nodes, not the number of edges
    } else if (G.isEmpty()) {
        vd = 0;
    } else {
        ConnectedComponents cc(G);
        DEBUG("finding connected components");
        cc.run();
        INFO("Number of components ", cc.numberOfComponents());
        DEBUG("Estimating size of the largest component");
        const auto compSizes = cc.getComponentSizes();
        vd = std::max_element(
                 compSizes.begin(), compSizes.end(),
                 [](const auto &c1, const auto &c2) -> bool { return c1.second < c2.second; })
                 ->second;

        INFO("Largest component size: ", vd);
    }
    return vd;
}

} /* namespace NetworKit */
