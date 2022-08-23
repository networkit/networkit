/*
 * dynSSSPGTest.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include <gtest/gtest.h>

#include <random>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/DynBFS.hpp>
#include <networkit/distance/DynDijkstra.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class DynSSSPGTest : public testing::Test {};

TEST_F(DynSSSPGTest, testDynamicBFS_1edge) {
    /* Graph:
        0    3   6
         \  / \ /
          2    5
         /  \ / \
        1    4   7
     */
    count n = 8;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    BFS bfs(G, 0);
    bfs.run();
    DynBFS dbfs(G, 0);
    dbfs.run();
    std::vector<GraphEvent> batch(1);
    batch[0].type = GraphEvent::EDGE_ADDITION;
    batch[0].u = 0;
    batch[0].v = 6;
    batch[0].w = 1.0;
    for (GraphEvent edge : batch) {
        G.addEdge(edge.u, edge.v, edge.w);
    }
    dbfs.updateBatch(batch);
    bfs.run();
    G.forNodes([&](node i) {
        EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
        EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
    });
}

TEST_F(DynSSSPGTest, testDynamicBFSEdgeDeletion) {
    /* Graph:
        0    3   6
         \  / \ /
          2    5
         /  \ / \
        1    4   7
     */
    Graph G(8);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    BFS bfs(G, 0);
    DynBFS dbfs(G, 0);
    dbfs.run();
    std::vector<GraphEvent> deletionBatch;
    deletionBatch.emplace_back(GraphEvent::EDGE_REMOVAL, 2, 3);
    for (GraphEvent edge : deletionBatch) {
        G.removeEdge(edge.u, edge.v);
    }
    dbfs.updateBatch(deletionBatch);
    bfs.run();
    G.forNodes([&](node i) {
        EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
        EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
    });
}

TEST_F(DynSSSPGTest, testDynamicBFS_batch) {
    /* Graph:
            0    3   6
            \  / \ /
                2    5
            /  \ / \
            1    4   7
    */
    count n = 8;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    BFS bfs(G, 0);
    bfs.run();
    DynBFS dbfs(G, 0);
    dbfs.run();
    std::vector<GraphEvent> batch(3);
    batch[0].type = GraphEvent::EDGE_ADDITION;
    batch[0].u = 3;
    batch[0].v = 7;
    batch[0].w = 1.0;
    batch[1].type = GraphEvent::EDGE_ADDITION;
    batch[1].u = 0;
    batch[1].v = 5;
    batch[1].w = 1.0;
    batch[2].type = GraphEvent::EDGE_ADDITION;
    batch[2].u = 2;
    batch[2].v = 7;
    batch[2].w = 1.0;
    for (GraphEvent edge : batch) {
        G.addEdge(edge.u, edge.v, edge.w);
    }
    dbfs.updateBatch(batch);
    bfs.run();
    G.forNodes([&](node i) {
        EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
        EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
    });
}

TEST_F(DynSSSPGTest, testDynamicBFS_batchDeletion) {
    /* Graph:
            0    3   6
            \  / \ /
              2 - 5
            /  \ / \
            1    4   7
    */
    Graph G(8);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);
    G.addEdge(2, 5);

    BFS bfs(G, 0);
    DynBFS dbfs(G, 0);
    dbfs.run();
    std::vector<GraphEvent> batch;
    batch.emplace_back(GraphEvent::EDGE_REMOVAL, 2, 3);
    batch.emplace_back(GraphEvent::EDGE_REMOVAL, 4, 5);
    batch.emplace_back(GraphEvent::EDGE_REMOVAL, 5, 7);
    for (GraphEvent edge : batch) {
        G.removeEdge(edge.u, edge.v);
    }
    dbfs.updateBatch(batch);
    bfs.run();
    G.forNodes([&](node i) {
        EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
        EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
    });
}

TEST_F(DynSSSPGTest, testDynamicDijkstra) {
    /* Graph:
       0    3   6
        \  / \ /
         2 -- 5
        /  \ / \
       1    4   7

       Edges in the upper row have weight 3,
       the edge in the middle row has weight 1.5,
       edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    Dijkstra dij(G, 0);
    dij.run();
    DynDijkstra ddij(G, 0);
    ddij.run();
    std::vector<GraphEvent> batch(3);
    batch[0].type = GraphEvent::EDGE_ADDITION;
    batch[0].u = 0;
    batch[0].v = 4;
    batch[0].w = 1.0;
    batch[1].type = GraphEvent::EDGE_ADDITION;
    batch[1].u = 1;
    batch[1].v = 4;
    batch[1].w = 1.0;
    batch[2].type = GraphEvent::EDGE_ADDITION;
    batch[2].u = 6;
    batch[2].v = 7;
    batch[2].w = 3.0;
    for (GraphEvent edge : batch) {
        G.addEdge(edge.u, edge.v, edge.w);
    }
    ddij.updateBatch(batch);
    dij.run();
    G.forNodes([&](node i) {
        EXPECT_EQ(dij.distance(i), ddij.distance(i));
        EXPECT_EQ(dij.numberOfPaths(i), ddij.numberOfPaths(i));
    });
}

TEST_F(DynSSSPGTest, testDynamicDijkstraDeletion) {
    /* Graph:
       0    3   6
        \  / \ /
         2 -- 5
        /  \ / \
       1    4   7

       Edges in the upper row have weight 3,
       the edge in the middle row has weight 1.5,
       edges in the lower row have weight 2.
    */
    Graph G(8, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    Dijkstra dij(G, 0);
    DynDijkstra ddij(G, 0);
    ddij.run();
    std::vector<GraphEvent> batch;
    batch.emplace_back(GraphEvent::EDGE_REMOVAL, 2, 3);
    batch.emplace_back(GraphEvent::EDGE_REMOVAL, 4, 5);
    batch.emplace_back(GraphEvent::EDGE_REMOVAL, 5, 7);
    for (GraphEvent edge : batch) {
        G.removeEdge(edge.u, edge.v);
    }
    ddij.updateBatch(batch);
    dij.run();
    G.forNodes([&](node i) {
        EXPECT_EQ(dij.distance(i), ddij.distance(i));
        EXPECT_EQ(dij.numberOfPaths(i), ddij.numberOfPaths(i));
    });
}

TEST_F(DynSSSPGTest, testDynamicBFSGeneratedGraph) {
    Aux::Random::setSeed(1, true);
    DorogovtsevMendesGenerator generator(500);
    Graph G = generator.generate();
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    DynBFS dyn_bfs(G, 0);
    BFS bfs(G, 0);
    dyn_bfs.run();
    bfs.run();
    DEBUG("Before the edge insertion: ");
    count nInsertions = 750, i = 0;
    while (i < nInsertions) {
        DEBUG("Sampling a new edge");
        node v1 = GraphTools::randomNode(G);
        node v2 = GraphTools::randomNode(G);
        if (v1 != v2 && !G.hasEdge(v1, v2)) {
            i++;
            DEBUG("Adding edge number ", i);
            G.addEdge(v1, v2);
            std::vector<GraphEvent> batch;
            batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
            DEBUG("Running update with dynamic bfs");
            dyn_bfs.updateBatch(batch);
            DEBUG("Running from scratch with bfs");
            bfs.run();
            G.forNodes([&](node i) {
                EXPECT_EQ(dyn_bfs.distance(i), bfs.distance(i));
                EXPECT_EQ(dyn_bfs.numberOfPaths(i), bfs.numberOfPaths(i));
            });
        }
    }
}

TEST_F(DynSSSPGTest, testDynamicBFSGeneratedGraphEdgeDeletion) {
    Aux::Random::setSeed(1, true);
    DorogovtsevMendesGenerator generator(500);
    Graph G = generator.generate();
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    DynBFS dyn_bfs(G, 0);
    BFS bfs(G, 0);
    dyn_bfs.run();
    DEBUG("Before the edge insertion: ");
    for (int i = 0; i < 750; ++i) {
        DEBUG("Sampling a new edge");
        auto randomEdge = GraphTools::randomEdge(G);
        DEBUG("Adding edge number ", i);
        G.removeEdge(randomEdge.first, randomEdge.second);
        std::vector<GraphEvent> batch;
        batch.emplace_back(GraphEvent::EDGE_REMOVAL, randomEdge.first, randomEdge.second);
        DEBUG("Running update with dynamic bfs");
        dyn_bfs.updateBatch(batch);
        DEBUG("Running from scratch with bfs");
        bfs.run();
        G.forNodes([&](node i) {
            EXPECT_EQ(dyn_bfs.distance(i), bfs.distance(i));
            EXPECT_EQ(dyn_bfs.numberOfPaths(i), bfs.numberOfPaths(i));
        });
    }
}

TEST_F(DynSSSPGTest, testDynamicDijkstraGeneratedGraph) {
    Aux::Random::setSeed(1, true);
    auto G = DorogovtsevMendesGenerator{1000}.generate();
    G = GraphTools::toWeighted(G);
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());

    DynDijkstra dyn_dij(G, 0);
    dyn_dij.run();

    Dijkstra dij(G, 0); // Baseline

    for (int i = 0; i < 10; ++i) {
        DEBUG("Sampling a new edge");
        node v1 = GraphTools::randomNode(G), v2 = GraphTools::randomNode(G);
        while (v1 == v2 || G.hasEdge(v1, v2))
            v1 = GraphTools::randomNode(G), v2 = GraphTools::randomNode(G);

        DEBUG("Adding edge number ", i);
        G.addEdge(v1, v2);
        DEBUG("Running update with dynamic dijkstra");
        dyn_dij.update(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, defaultEdgeWeight));
        DEBUG("Running from scratch with dijkstra");
        dij.run();
        G.forNodes([&](node i) {
            EXPECT_DOUBLE_EQ(dyn_dij.distance(i), dij.distance(i));
            EXPECT_EQ(dyn_dij.numberOfPaths(i), dij.numberOfPaths(i));
        });
    }
}

TEST_F(DynSSSPGTest, testDynamicDijkstraGeneratedGraphEdgeDeletion) {
    Aux::Random::setSeed(1, true);
    DorogovtsevMendesGenerator generator(1000);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    DynDijkstra dyn_dij(G, 0);
    Dijkstra dij(G, 0);
    dyn_dij.run();
    DEBUG("Before the edge deletion: ");
    for (count i = 0; i < 10; i++) {
        DEBUG("Selecting a random edge");
        auto randomEdge = GraphTools::randomEdge(G);
        DEBUG("Deleting edge number ", i);
        G.removeEdge(randomEdge.first, randomEdge.second);
        std::vector<GraphEvent> batch;
        batch.emplace_back(GraphEvent::EDGE_REMOVAL, randomEdge.first, randomEdge.second);
        DEBUG("Running update with dynamic dijkstra");
        dyn_dij.updateBatch(batch);
        DEBUG("Running from scratch with dijkstra");
        dij.run();
        G.forNodes([&](node i) {
            EXPECT_EQ(dyn_dij.distance(i), dij.distance(i));
            EXPECT_EQ(dyn_dij.numberOfPaths(i), dij.numberOfPaths(i));
        });
    }
}

TEST_F(DynSSSPGTest, testDynamicDijkstraBatches) {
    Aux::Random::setSeed(1, true);
    std::default_random_engine random_generator;
    std::normal_distribution<double> distribution(100, 10);
    DorogovtsevMendesGenerator generator(100);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    // add random normal weights to G

    G.forNodes([&](node source) {
        DynDijkstra dyn_dij(G, source, true);
        Dijkstra dij(G, source);
        dyn_dij.run();
        dij.run();
        DEBUG("Before the edge insertion: ");
        count batchSize = 8;
        count nBatches = 1, i = 0;
        for (count j = 0; j < nBatches; j++) {
            std::vector<GraphEvent> batch;
            i = 0;
            while (i < batchSize) {
                DEBUG("Sampling a new edge");
                node v1 = GraphTools::randomNode(G);
                node v2 = GraphTools::randomNode(G);
                if (v1 != v2 && !G.hasEdge(v1, v2)) {
                    i++;
                    double number = distribution(random_generator);
                    G.addEdge(v1, v2, number);
                    batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, number));
                }
            }
            DEBUG("batch size: ", batch.size());
            DEBUG("Updating with dynamic dijkstra");
            dyn_dij.updateBatch(batch);
            DEBUG("Running from scratch with dijkstra");
            dij.run();
            G.forNodes([&](node i) {
                EXPECT_EQ(dyn_dij.distance(i), dij.distance(i));
                EXPECT_EQ(dyn_dij.numberOfPaths(i), dij.numberOfPaths(i));
                if (i != source) {
                    ASSERT_NE(dyn_dij.distance(i), 0);
                }
            });
        }
    });
}

TEST_F(DynSSSPGTest, testDynamicDijkstraBatchesEdgeDeletions) {
    Aux::Random::setSeed(1, true);
    ErdosRenyiGenerator generator(100, 0.25, false, false);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    // delete random edges from G

    G.forNodes([&](node source) {
        DynDijkstra dyn_dij(G, source, true);
        Dijkstra dij(G, source);
        dyn_dij.run();
        DEBUG("Before the edge deletion: ");
        std::vector<GraphEvent> batch;
        for (count i = 0; i < 8; ++i) {
            DEBUG("Selecting a random edge");
            auto randomEdge = GraphTools::randomEdge(G);
            G.removeEdge(randomEdge.first, randomEdge.second);
            batch.emplace_back(GraphEvent::EDGE_REMOVAL, randomEdge.first, randomEdge.second);
        }
        DEBUG("batch size: ", batch.size());
        DEBUG("Updating with dynamic dijkstra");
        dyn_dij.updateBatch(batch);
        DEBUG("Running from scratch with dijkstra");
        dij.run();
        G.forNodes([&](node i) {
            EXPECT_EQ(dyn_dij.distance(i), dij.distance(i));
            EXPECT_EQ(dyn_dij.numberOfPaths(i), dij.numberOfPaths(i));
            if (i != source) {
                EXPECT_NE(dyn_dij.distance(i), 0);
            }
        });
    });
}

} /* namespace NetworKit */
