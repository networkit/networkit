/*
 * DistanceGTest.cpp
 *
 *  Created on: Sep 04, 2015
 *      Author: Maximilian Vogel
 */

#include <limits>
#include <gtest/gtest.h>

#include <networkit/distance/APSP.hpp>
#include <networkit/distance/AStar.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/BidirectionalBFS.hpp>
#include <networkit/distance/BidirectionalDijkstra.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/EffectiveDiameter.hpp>
#include <networkit/distance/EffectiveDiameterApproximation.hpp>
#include <networkit/distance/HopPlotApproximation.hpp>
#include <networkit/distance/IncompleteDijkstra.hpp>
#include <networkit/distance/MultiTargetBFS.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/distance/NeighborhoodFunction.hpp>
#include <networkit/distance/NeighborhoodFunctionApproximation.hpp>
#include <networkit/distance/NeighborhoodFunctionHeuristic.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/distance/SPSP.hpp>

#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {
class DistanceGTest : public testing::TestWithParam<std::pair<bool, bool>> {

protected:
    static constexpr auto infdist = std::numeric_limits<edgeweight>::max();

    bool isDirected() const noexcept;
    bool isWeighted() const noexcept;

    Graph generateERGraph(count n, double p) const {
        auto G = ErdosRenyiGenerator(n, p, isDirected()).generate();
        if (isWeighted()) {
            G = GraphTools::toWeighted(G);
            G.forEdges([&G](node u, node v) { G.setWeight(u, v, Aux::Random::probability()); });
        }

        return G;
    }
};

constexpr edgeweight DistanceGTest::infdist;

INSTANTIATE_TEST_SUITE_P(InstantiationName, DistanceGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

bool DistanceGTest::isWeighted() const noexcept {
    return GetParam().first;
}

bool DistanceGTest::isDirected() const noexcept {
    return GetParam().second;
}

TEST_F(DistanceGTest, testVertexDiameterPedantically) {
    DorogovtsevMendesGenerator generator(1000);
    Graph G1 = generator.generate();
    Graph G = Graph(G1, true, false);
    Diameter diam(G, DiameterAlgo::ESTIMATED_PEDANTIC);
    diam.run();
    count vd = diam.getDiameter().first;
    EXPECT_EQ(1000, vd);
}

TEST_F(DistanceGTest, testAStar) {

    // Builds a mesh graph with the given number of rows and columns
    auto buildMesh = [](count rows, count cols) -> Graph {
        Graph G(rows * cols, false, false);

        for (count i = 0; i < rows; ++i) {
            for (count j = 0; j < cols; ++j) {
                if (j < cols - 1) {
                    G.addEdge(i * cols + j, i * cols + j + 1);
                }
                if (i < rows - 1) {
                    G.addEdge(i * cols + j, (i + 1) * cols + j);
                }
            }
        }

        return G;
    };

    // Test AStar on a mesh for every pair of nodes (u, v), u != v.
    auto testMesh = [&](count rows, count cols) {
        auto G = buildMesh(rows, cols);

        // Test AStar for every pair.
        G.forNodePairs([&](node source, node target) {
            // Three distance heuristics
            std::vector<double> zeroDist(rows * cols, 0);
            std::vector<double> exactDist(rows * cols);
            std::vector<double> eucledianDist(rows * cols);
            for (node u = 0; u < rows * cols; ++u) {
                count rowU = u / cols;
                count colU = u % cols;
                count rowT = target / cols;
                count colT = target % cols;
                exactDist[u] = static_cast<double>((rowU > rowT ? rowU - rowT : rowT - rowU)
                                                   + (colU > colT ? colU - colT : colT - colU));
                double rowDiff =
                    std::abs(static_cast<double>(u / cols) - static_cast<double>(rowT));
                double colDiff =
                    std::abs(static_cast<double>(u % cols) - static_cast<double>(rowT));
                eucledianDist[u] = std::sqrt(std::pow(rowDiff, 2) + std::pow(colDiff, 2));
            };
            BFS bfs(G, source, true, false, target);
            bfs.run();

            for (const auto &heu : {zeroDist, exactDist, eucledianDist}) {
                AStar astar(G, heu, source, target, true);
                astar.run();
                EXPECT_EQ(astar.getDistance(), bfs.distance(target));

                auto path = astar.getPath();
                EXPECT_EQ(path.size(), bfs.getPath(target).size() - 2);

                if (!path.size()) {
                    return;
                }
                for (size_t i = 0; i < path.size() - 1; ++i) {
                    EXPECT_TRUE(G.hasEdge(path[i], path[i + 1]));
                }
            }
        });
    };

    // Test cases
    testMesh(10, 10);
    testMesh(1, 15);
    testMesh(25, 5);
}

TEST_P(DistanceGTest, testIncompleteDijkstra) {
    Aux::Random::setSeed(42, false);
    auto G = generateERGraph(500, 0.05);

    G.forNodes([&](const node source) {
        Dijkstra dij(G, source, false, false);
        dij.run();
        const auto dists = dij.getDistances();
        const count reachable = std::count_if(
            dists.begin(), dists.end(), [](const edgeweight dist) { return dist != infdist; });

        const std::vector<node> sources({source});
        IncompleteDijkstra iDij(&G, sources);

        for (count i = 0; i < reachable; ++i) {
            EXPECT_TRUE(iDij.hasNext());
            const auto next = iDij.next();
            EXPECT_DOUBLE_EQ(next.second, dists[next.first]);
        }

        EXPECT_FALSE(iDij.hasNext());
    });
}

TEST_P(DistanceGTest, testBidirectionalBFS) {
    Aux::Random::setSeed(42, false);
    auto G = generateERGraph(500, isDirected() ? 0.05 : 0.02);

    for (int i = 0; i < 10; ++i) {
        node source = GraphTools::randomNode(G);
        node target = GraphTools::randomNode(G);
        while (source == target)
            target = GraphTools::randomNode(G);
        BFS bfs(G, source, true, false, target);
        bfs.run();
        BidirectionalBFS bbfs(G, source, target, true);
        bbfs.run();

        if (bfs.distance(target) < G.upperNodeIdBound()) {
            // Source reaches target
            EXPECT_EQ(bfs.distance(target), bbfs.getDistance());
            const auto path = bbfs.getPath();
            EXPECT_EQ(path.size(), bbfs.getDistance() - 1);

            if (!path.empty()) { // At least two edges in the path
                EXPECT_TRUE(G.hasEdge(source, path.front()));
                EXPECT_TRUE(G.hasEdge(path.back(), target));
                for (count i = 1; i < path.size() - 1; ++i)
                    EXPECT_TRUE(G.hasEdge(path[i], path[i + 1]));
            } else
                EXPECT_TRUE(G.hasEdge(source, target));
        } else
            EXPECT_EQ(infdist, bbfs.getDistance());
    };
}

TEST_P(DistanceGTest, testBidirectionalDijkstra) {
    Aux::Random::setSeed(42, false);
    auto G = generateERGraph(500, isDirected() ? 0.05 : 0.02);

    for (int i = 0; i < 10; ++i) {
        node source = GraphTools::randomNode(G);
        node target = GraphTools::randomNode(G);
        while (source == target)
            target = GraphTools::randomNode(G);
        BidirectionalDijkstra bdij(G, source, target, true);
        bdij.run();
        Dijkstra dij(G, source, true, true, target);
        dij.run();

        EXPECT_DOUBLE_EQ(bdij.getDistance(), dij.distance(target));
        const auto path = bdij.getPath();

        if (dij.distance(target) != infdist) {
            // Source reaches target
            EXPECT_DOUBLE_EQ(dij.distance(target), bdij.getDistance());

            if (!path.empty()) { // At least two edges in the path
                EXPECT_TRUE(G.hasEdge(source, path.front()));
                for (size_t i = 0; i < path.size() - 1; ++i)
                    EXPECT_TRUE(G.hasEdge(path[i], path[i + 1]));
                EXPECT_TRUE(G.hasEdge(path.back(), target));
            } else
                EXPECT_TRUE(G.hasEdge(source, target));
        } else
            EXPECT_EQ(infdist, bdij.getDistance());
    };
}

TEST_F(DistanceGTest, testExactDiameter) {
    using namespace std;

    vector<pair<string, count>> testInstances = {pair<string, count>("lesmis", 14),
                                                 pair<string, count>("jazz", 6),
                                                 pair<string, count>("celegans_metabolic", 7)};

    for (auto testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance.first + ".graph");
        Diameter diam(G, DiameterAlgo::EXACT);
        diam.run();
        count diameter = diam.getDiameter().first;
        EXPECT_EQ(diameter, testInstance.second);
    }
}

TEST_F(DistanceGTest, testEstimatedDiameterRange) {
    using namespace std;

    vector<pair<string, count>> testInstances = {pair<string, count>("celegans_metabolic", 7),
                                                 pair<string, count>("jazz", 6)};

    for (auto testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance.first + ".graph");
        Diameter diam(G, DiameterAlgo::ESTIMATED_RANGE, 0.1);
        diam.run();
        std::pair<count, count> range = diam.getDiameter();
        EXPECT_GE(testInstance.second, range.first);
        EXPECT_LE(testInstance.second, range.second);
    }
}
TEST_F(DistanceGTest, testPedanticDiameterErdos) {
    count n = 5000;
    ErdosRenyiGenerator gen(n, 0.001);
    Graph G1 = gen.generate();
    Diameter diam(G1, DiameterAlgo::ESTIMATED_PEDANTIC);
    diam.run();
    count diameter = diam.getDiameter().first;
    ASSERT_LE(diameter, n);
}

TEST_F(DistanceGTest, testEffectiveDiameterMinimal) {
    // Minimal example from the paper
    Graph G(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(4, 0);
    EffectiveDiameterApproximation aef(G);
    aef.run();
    double effective = aef.getEffectiveDiameter();
    Diameter diam(G, DiameterAlgo::EXACT);
    diam.run();
    count exact = diam.getDiameter().first;
    EXPECT_LE(effective, exact);
}

TEST_F(DistanceGTest, testEffectiveDiameter) {

    using namespace std;

    vector<string> testInstances = {"celegans_metabolic", "jazz", "lesmis"};

    for (auto testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance + ".graph");
        EffectiveDiameterApproximation aef(G);
        aef.run();
        double effective = aef.getEffectiveDiameter();
        Diameter diam(G, DiameterAlgo::EXACT);
        diam.run();
        count exact = diam.getDiameter().first;
        EXPECT_LE(effective, exact);
    }
}

TEST_F(DistanceGTest, testEffectiveDiameterExact) {

    using namespace std;

    vector<string> testInstances = {"celegans_metabolic", "jazz", "lesmis"};

    for (auto testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance + ".graph");
        EffectiveDiameter ed(G);
        ed.run();
        double effective = ed.getEffectiveDiameter();
        Diameter diam(G, DiameterAlgo::EXACT);
        diam.run();
        count exact = diam.getDiameter().first;
        EXPECT_LE(effective, exact);
    }

    const double tol = 1e-3;

    /* Graph: n=20, threshold: 20*0.9 = 18 nodes
        1--3--5--7---9
        |  |  |  |   |
        2--4--6--8--10
            |     |
            11----12
                |
            13--14--15
                |
            18--16--17--19
                    |
                    20
    Number of steps needed per node: (1-20)
    (7+6+6+5+6+5+5+4+6+5+4+4+5+4+5+5+6+6+7+7) / 20 = 5.4
    */
    count n1 = 20;
    Graph G1(n1);

    G1.addEdge(0, 1);
    G1.addEdge(0, 2);
    G1.addEdge(1, 3);
    G1.addEdge(2, 3);
    G1.addEdge(2, 4);
    G1.addEdge(3, 5);
    G1.addEdge(3, 10);
    G1.addEdge(4, 5);
    G1.addEdge(4, 6);
    G1.addEdge(5, 7);
    G1.addEdge(6, 8);
    G1.addEdge(6, 7);
    G1.addEdge(7, 9);
    G1.addEdge(7, 11);
    G1.addEdge(8, 9);
    G1.addEdge(10, 11);
    G1.addEdge(11, 13);
    G1.addEdge(12, 13);
    G1.addEdge(13, 14);
    G1.addEdge(13, 15);
    G1.addEdge(15, 16);
    G1.addEdge(15, 17);
    G1.addEdge(16, 18);
    G1.addEdge(16, 19);

    EffectiveDiameter ed(G1);
    ed.run();
    double effective1 = ed.getEffectiveDiameter();
    EXPECT_NEAR(5.4, effective1, tol);

    /* Graph: n=21, threshold: 21*0.9 = 18.9 => 19 nodes
                13---------------3
                    |               |
                ---14--12--|        |
                |   |   |  |        |
    1--21--18--16--15   |  |        |
                |       |  |        |
        20--17------10--8        |
                |       |  |        |
            19       9--7--5--6--4--11
                                    |
                                    2
Number of steps needed per node: (1-21)
(8+7+5+6+6+6+5+5+5+5+7+5+4+4+5+5+5+6+6+6+7) / 21 = 5.619047
*/
    count n2 = 21;
    Graph G2(n2);

    G2.addEdge(0, 20);
    G2.addEdge(1, 3);
    G2.addEdge(2, 3);
    G2.addEdge(2, 12);
    G2.addEdge(3, 5);
    G2.addEdge(3, 10);
    G2.addEdge(4, 5);
    G2.addEdge(4, 6);
    G2.addEdge(6, 7);
    G2.addEdge(6, 8);
    G2.addEdge(7, 9);
    G2.addEdge(7, 11);
    G2.addEdge(8, 9);
    G2.addEdge(9, 11);
    G2.addEdge(9, 16);
    G2.addEdge(11, 13);
    G2.addEdge(12, 13);
    G2.addEdge(13, 14);
    G2.addEdge(13, 15);
    G2.addEdge(14, 15);
    G2.addEdge(15, 16);
    G2.addEdge(15, 17);
    G2.addEdge(16, 18);
    G2.addEdge(16, 19);
    G2.addEdge(17, 20);

    EffectiveDiameter ed2(G2);
    ed2.run();
    double effective2 = ed2.getEffectiveDiameter();
    EXPECT_NEAR(5.619047, effective2, tol);
}

TEST_F(DistanceGTest, testHopPlotApproximation) {
    using namespace std;

    vector<string> testInstances = {"celegans_metabolic", "lesmis"};

    const double tol = 1e-2;

    for (auto &testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance + ".graph");
        HopPlotApproximation hp(G);
        hp.run();
        map<count, double> hopPlot = hp.getHopPlot();
        for (count i = 1; i < hopPlot.size(); i++) {
            EXPECT_LE(hopPlot[i - 1], hopPlot[i] + tol);
        }
    }
}

TEST_F(DistanceGTest, testNeighborhoodFunctionApproximation) {
    METISGraphReader reader;
    Graph G = GraphTools::toUnweighted(reader.read("input/lesmis.graph"));
    NeighborhoodFunction nf(G);
    nf.run();
    auto exact = nf.getNeighborhoodFunction();
    NeighborhoodFunctionApproximation anf(G);
    anf.run();
    auto approximated = anf.getNeighborhoodFunction();
    EXPECT_EQ(exact.size(), approximated.size());
}

TEST_F(DistanceGTest, testNeighborhoodFunctionHeuristic) {
    METISGraphReader reader;
    Graph G = GraphTools::toUnweighted(reader.read("input/lesmis.graph"));
    NeighborhoodFunction nf(G);
    nf.run();
    auto exact = nf.getNeighborhoodFunction();
    NeighborhoodFunctionHeuristic anf(G);
    anf.run();
    auto heuristic = anf.getNeighborhoodFunction();
    EXPECT_EQ(exact.size(), heuristic.size());
}

TEST_P(DistanceGTest, testSPSP) {
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = generateERGraph(100, 0.15);

        APSP apsp(G);
        apsp.run();
        const auto gt = apsp.getDistances();

        for (count nSources : {1, 10, 50, 100}) {
            auto sources = GraphTools::randomNodes(G, nSources);

            SPSP spsp(G, sources.begin(), sources.end());
            spsp.run();
            const auto nodemap = spsp.getSourceIndexMap();
            const auto distances = spsp.getDistances();

            for (node source : sources) {
                const auto &curDist = distances[nodemap.at(source)];
                G.forNodes([&](node target) {
                    EXPECT_DOUBLE_EQ(gt[source][target], spsp.getDistance(source, target));
                    EXPECT_DOUBLE_EQ(gt[source][target], curDist[target]);
                });
            }
        }
    }
}

TEST_P(DistanceGTest, testSPSPWithTargets) {
    Aux::Random::setSeed(42, true);
    const auto G = generateERGraph(100, 0.15);

    APSP apsp(G);
    apsp.run();

    for (count nSources : {1, 10, 100}) {
        const auto sources = GraphTools::randomNodes(G, nSources);
        for (count nTargets : {1, 50, 100}) {
            const auto targets = GraphTools::randomNodes(G, nTargets);

            SPSP spsp(G, sources.begin(), sources.end(), targets.begin(), targets.end());
            spsp.run();
            const auto sourceMap = spsp.getSourceIndexMap(), targetMap = spsp.getTargetIndexMap();
            const auto distances = spsp.getDistances();
            EXPECT_EQ(distances.size(), nSources);
            EXPECT_EQ(distances.front().size(), nTargets);

            for (node source : sources)
                for (node target : targets)
                    EXPECT_EQ(apsp.getDistance(source, target), spsp.getDistance(source, target));
        }
    }
}

TEST_P(DistanceGTest, testSPSPWithUnreachableTarget) {
    Aux::Random::setSeed(42, true);
    auto G = generateERGraph(100, 0.15);
    const node unreachable = G.addNode();

    APSP apsp(G);
    apsp.run();

    const std::vector<node> sources = {0}, targets = {unreachable};
    SPSP spsp(G, sources.begin(), sources.end(), targets.begin(), targets.end());
    spsp.run();

    for (node source : sources)
        for (node target : targets)
            EXPECT_DOUBLE_EQ(apsp.getDistance(source, target), spsp.getDistance(source, target));
}

TEST_P(DistanceGTest, testMultiTargetBFS) {
    Aux::Random::setSeed(42, true);
    const auto G = generateERGraph(100, 0.15);
    const auto source = GraphTools::randomNode(G);

    BFS bfs(G, source, false);
    bfs.run();
    const auto distances = bfs.getDistances();

    for (count nTargets : {1, 10, 20}) {
        const auto targets = GraphTools::randomNodes(G, nTargets);
        MultiTargetBFS mtBFS(G, source, targets.begin(), targets.end());
        mtBFS.run();
        const auto tgtDists = mtBFS.getDistances();
        const auto tgtIdx = mtBFS.getTargetIndexMap();

        for (node target : targets)
            EXPECT_EQ(distances[target], tgtDists[tgtIdx.at(target)]);
    }
}

TEST_P(DistanceGTest, testMultiTargetDijkstra) {
    Aux::Random::setSeed(42, true);
    const auto G = generateERGraph(100, 0.15);
    const auto source = GraphTools::randomNode(G);

    Dijkstra dij(G, source, false);
    dij.run();
    const auto distances = dij.getDistances();

    for (count nTargets : {1, 10, 20}) {
        const auto targets = GraphTools::randomNodes(G, nTargets);
        MultiTargetDijkstra mtDij(G, source, targets.begin(), targets.end());
        mtDij.run();
        const auto tgtDists = mtDij.getDistances();
        const auto tgtIdx = mtDij.getTargetIndexMap();

        for (node target : targets)
            EXPECT_DOUBLE_EQ(distances[target], tgtDists[tgtIdx.at(target)]);
    }
}

TEST_P(DistanceGTest, testPrunedLandmarkLabeling) {
    Aux::Random::setSeed(42, false);
    Graph G = ErdosRenyiGenerator{500, 0.01, isDirected()}.generate();
    PrunedLandmarkLabeling pll(G);
    pll.run();

    ASSERT_TRUE(pll.hasFinished());

    APSP apsp(G);
    apsp.run();

    const count infDistInt = std::numeric_limits<count>::max();
    const double infDistDouble = std::numeric_limits<double>::max();

    G.forNodePairs([&](node u, node v) {
        double distUV = apsp.getDistance(u, v);
        if (distUV == infDistDouble)
            EXPECT_EQ(pll.query(u, v), infDistInt);
        else
            EXPECT_EQ(pll.query(u, v), distUV);
    });
}

} /* namespace NetworKit */
