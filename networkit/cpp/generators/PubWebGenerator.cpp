/*
 * PubWebGenerator.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */
// networkit-format

#include <cmath>
#include <queue>
#include <set>

#include <networkit/auxiliary/Random.hpp>

#include <networkit/generators/PubWebGenerator.hpp>

namespace NetworKit {

PubWebGenerator::PubWebGenerator(count n, count numDenseAreas, coord neighRad,
                                 count maxNumNeighbors)
    : n(n), numDenseAreas(numDenseAreas), neighRad(neighRad), maxNeigh(maxNumNeighbors) {}

coord2d PubWebGenerator::intoUnitSquare(coord2d pt) const noexcept {
    auto toUnitSquare = [](coord z) {
        if (z > 1.0)
            return z - 1.0;
        if (z < 0.0)
            return z + 1.0;
        return z;
    };

    return {toUnitSquare(pt.first), toUnitSquare(pt.second)};
}

coord PubWebGenerator::squaredDistanceInUnitTorus(coord2d pt1, coord2d pt2) const noexcept {
    auto adjustForUnitTorus = [](coord z) -> coord {
        if (z > 0.5)
            return 1.0 - z;
        if (z < -0.5)
            return z + 1.0;
        return z;
    };

    coord distx = adjustForUnitTorus(pt1.first - pt2.first);
    coord disty = adjustForUnitTorus(pt1.second - pt2.second);

    return distx * distx + disty * disty;
}

// TODO: use ANN or similar library with appropriate space-partitioning data structure to
//       get rid of quadratic time complexity
void PubWebGenerator::determineNeighbors(Graph &g) {
    coord sqrNeighRad = neighRad * neighRad;

    using edge = std::pair<node, node>;
    std::set<std::pair<node, node>> eligibleEdges;

    auto isInRange([&](coord squaredDistance) { return (squaredDistance <= sqrNeighRad); });

    g.forNodes([&](node u) {
        std::priority_queue<std::pair<coord, edge>> pq;
        const auto p1 = coordinates[u];

        // fill PQ with neighbors in range
        g.forNodes([&](node v) {
            const auto sqrDist = squaredDistanceInUnitTorus(p1, coordinates[v]);

            if (isInRange(sqrDist)) {
                edge e = std::make_pair(std::min(u, v), std::max(u, v));
                pq.push(std::make_pair(-sqrDist, e));
            }
        });

        // mark maxNeigh nearest neighbors as eligible or insert them into graph g
        count end = std::min(maxNeigh, (count)pq.size());
        for (index i = 0; i < end; ++i) {
            const auto currentBest = pq.top();
            pq.pop();

            if (eligibleEdges.count(currentBest.second) > 0) {
                // edge is already marked => insert it
                edgeweight ew = BASE_WEIGHT / -currentBest.first;
                g.addEdge(currentBest.second.first, currentBest.second.second, ew);
                //				TRACE("add edge " , currentBest.second.first , "/" ,
                // currentBest.second.second);
            } else {
                // mark edge as eligible
                eligibleEdges.insert(currentBest.second);
            }
        }
    });
}

void PubWebGenerator::addNodesToArea(index area, count num, Graph &g) {
    for (index j = 0; j < num; ++j) {
        // compute random angle between [0, 2pi) and distance between [0, width/2]
        coord angle = Aux::Random::real() * 2.0 * PI;
        coord dist = Aux::Random::real() * denseAreaXYR[area].rad;

        // compute coordinates and adjust them
        coord x = denseAreaXYR[area].x + std::cos(angle) * dist;
        coord y = denseAreaXYR[area].y + std::sin(angle) * dist;

        // create vertex with these coordinates
        g.addNode();
        coordinates.push_back(intoUnitSquare({x, y}));
    }
}

void PubWebGenerator::chooseDenseAreaSizes() {
    denseAreaXYR.resize(numDenseAreas);

    for (index area = 0; area < numDenseAreas; ++area) {
        // anti-quadratic probability distribution
        coord f = Aux::Random::real() * MIN_MAX_DENSE_AREA_FACTOR + 1.0;
        denseAreaXYR[area].rad = (MAX_DENSE_AREA_RADIUS * f * f)
                                 / (MIN_MAX_DENSE_AREA_FACTOR * MIN_MAX_DENSE_AREA_FACTOR);
    }
}

// compute number of nodes per cluster, each cluster has approx. same density
void PubWebGenerator::chooseClusterSizes() {
    auto f = std::accumulate(denseAreaXYR.begin(), denseAreaXYR.end(), 0.0,
                             [](coord sum, circle c) { return sum + pow(c.rad, 1.5); });

    f = (n * (numDenseAreas / (numDenseAreas + 2.0))) / f;

    numPerArea.reserve(numDenseAreas);
    for (auto circle : denseAreaXYR) {
        numPerArea.emplace_back(std::round(f * pow(circle.rad, 1.5)));
    }
}

void PubWebGenerator::fillDenseAreas(Graph &g) {
    for (index area = 0; area < numDenseAreas; ++area) {
        // choose center randomly, ensure complete cluster is within (0,1) without modifications
        denseAreaXYR[area].x = Aux::Random::real();
        denseAreaXYR[area].y = Aux::Random::real();
        addNodesToArea(area, numPerArea[area], g);
    }
}

// randomly spread remaining vertices over whole area
void PubWebGenerator::spreadRemainingNodes(Graph &g) {
    while (g.numberOfNodes() < n) {
        g.addNode(); // TODO: Replace with g.addNodes()
        coordinates.emplace_back(Aux::Random::real(), Aux::Random::real());
    }
}

Graph PubWebGenerator::generate() {
    // init
    Graph G(0, true);

    // add vertices according to PubWeb distribution
    chooseDenseAreaSizes();
    chooseClusterSizes();
    fillDenseAreas(G);
    spreadRemainingNodes(G);
    determineNeighbors(G);

    G.shrinkToFit();
    return G;
}

} /* namespace NetworKit */
