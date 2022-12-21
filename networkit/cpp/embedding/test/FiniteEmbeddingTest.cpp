/*
 * FiniteEmbeddingTest.cpp
 *
 *  Created on: 16.10.2020
 *      Author: Klaus Ahrens
 */

#include <iomanip>
#include <iostream>

#include <gtest/gtest.h>

#include <networkit/embedding/Node2Vec.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class FiniteEmbeddingTest : public testing::Test {};

using Embeddings = std::vector<std::vector<float>>;

bool allFinite(const Embeddings &e) {
    for (auto &emb : e) {
        for (auto f : emb) {
            if (std::isinf(f)) {
                return false;
            }
        }
    }
    return true;
}

TEST_F(FiniteEmbeddingTest, testNode2VecOnDirectedGraph) {
    // Directed graph with two nodes.
    Graph G(2, false, true);
    G.addEdge(0, 1);
    EXPECT_THROW(NetworKit::Node2Vec{G}, std::runtime_error);
}

TEST_F(FiniteEmbeddingTest, testNode2VecOnGraphWithIsolatedNodes) {
    // Undirected graph with three nodes; node with id 2 is isolated.
    Graph G(3);
    G.addEdge(0, 1);
    EXPECT_THROW(NetworKit::Node2Vec{G}, std::runtime_error);
}

TEST_F(FiniteEmbeddingTest, testNode2VecOnGraphWithNonContinuousIds) {
    // Undirected graph with two nodes; only id 0 and id 2 are present.
    Graph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.removeNode(1);
    EXPECT_THROW(NetworKit::Node2Vec{G}, std::runtime_error);
}

TEST_F(FiniteEmbeddingTest, testFiniteEmbeddingOnSmallGraph) {

    NetworKit::METISGraphReader reader;

    NetworKit::Graph graph = reader.read("input/karate.graph");

    auto algo = NetworKit::Node2Vec(graph, .3, 1.3);

    algo.run();

    auto features = algo.getFeatures();

    EXPECT_TRUE(allFinite(features));
}

} // namespace NetworKit
