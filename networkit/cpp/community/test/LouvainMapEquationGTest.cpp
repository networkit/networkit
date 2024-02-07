/*
 * LouvainMapEquationGTest.cpp
 *
 * Created on: 2019-10-30
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class MapEquationGTest : public testing::Test {
public:
    void SetUp() { Aux::Random::setSeed(435913, false); }
};

void addClique(Graph &graph, Partition &groundTruth, node lowestId, node highestId) {
    index subsetId = groundTruth.upperBound();
    groundTruth.setUpperBound(subsetId + 1);
    for (node i = lowestId; i <= highestId; ++i) {
        groundTruth.addToSubset(subsetId, i);
        for (node j = i + 1; j <= highestId; ++j) {
            graph.addEdge(i, j);
        }
    }
}

TEST_F(MapEquationGTest, testLocalMoveSmall) {
    Aux::Random::setSeed(2342556, false);
    Graph G(10);
    Partition groundTruth(10);
    addClique(G, groundTruth, 0, 4);
    addClique(G, groundTruth, 5, 9);
    G.addEdge(0, 9);
    G.addEdge(1, 8);
    G.addEdge(2, 7);

    LouvainMapEquation mapequation(G, false, 256, "none");
    mapequation.run();
    auto partition = mapequation.getPartition();

    EXPECT_EQ(partition.getSubsets(), groundTruth.getSubsets());
}

TEST_F(MapEquationGTest, testLocalMove) {
    Aux::Random::setSeed(2342556, false);

    count n = 300;
    BarabasiAlbertGenerator generator(50, n);
    Graph G = generator.generate();

    // Generates an std::unordered_set<node> with the sequence of consecutive integers in [start,
    // end).
    auto generateConsecutiveNodes = [](node start, node end) -> std::unordered_set<node> {
        std::vector<node> nodes(end - start);
        std::iota(nodes.begin(), nodes.end(), start);
        return {nodes.begin(), nodes.end()};
    };

    Graph G1 = GraphTools::subgraphFromNodes(G, generateConsecutiveNodes(0, 100));
    const Graph G2 = GraphTools::subgraphFromNodes(G, generateConsecutiveNodes(100, 200));
    const Graph G3 = GraphTools::subgraphFromNodes(G, generateConsecutiveNodes(200, 300));

    // generate 3 independent clusterings (=ground truth)
    Partition groundTruth(n);
    groundTruth.setUpperBound(3);
    for (node i = 0; i < n; ++i) {
        if (i < 100) {
            groundTruth.addToSubset(0, i);
        } else if (100 <= i && i < 200) {
            groundTruth.addToSubset(1, i);
        } else {
            groundTruth.addToSubset(2, i);
        }
    }

    // G1 becomes clustered version of G
    GraphTools::merge(G1, G2);
    GraphTools::merge(G1, G3);

    LouvainMapEquation algo(G1, true, 32, "synchronous");
    algo.run();
    auto partition = algo.getPartition();

    EXPECT_EQ(partition.getSubsets(), groundTruth.getSubsets());
}

} // namespace NetworKit
