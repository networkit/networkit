/*
 * AttributeTest.cpp
 *
 *  Created on: 26.02.2022
 *      Author: Klaus Ahrens
 */

#include <cstdio>
#include <iostream>
#include <numeric>

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class AttributeTest : public testing::Test {};

TEST_F(AttributeTest, testSetGetOnExistingNodes) {

    Graph graph(15);

    const std::string name = "some int attribute";
    auto intAttr = graph.nodeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.nodeAttributes().find(name)->second->getName(), name);

    // set / get
    graph.forNodes([&](node n) { intAttr.set(n, int(n)); });
    EXPECT_EQ(intAttr.size(), graph.numberOfNodes());

    graph.forNodes([&](node n) { EXPECT_EQ(intAttr.get(n), int(n)); });
    EXPECT_EQ(intAttr.size(), graph.numberOfNodes());

    // by index proxy
    graph.forNodes([&](node n) { intAttr[n] = int(n); });
    EXPECT_EQ(intAttr.size(), graph.numberOfNodes());

    graph.forNodes([&](node n) { EXPECT_EQ(intAttr[n], int(n)); });
    EXPECT_EQ(intAttr.size(), graph.numberOfNodes());

    // Test iterator
    node u = 0;
    int att = 0;
    for (auto pair : intAttr) {
        EXPECT_EQ(pair.first, u);
        EXPECT_EQ(pair.second, att);
        ++u, ++att;
    }
}

TEST_F(AttributeTest, testReadWrite) {
    const int n = 10;
    const std::string path = "output/attribute.txt";
    const std::string name = "attribute";

    std::vector<double> values(n);
    std::iota(values.begin(), values.end(), 0);

    {
        Graph graph(n);
        auto attr = graph.nodeAttributes().attach<double>(name);
        graph.forNodes([&](node u) { attr.set(u, static_cast<double>(values[u])); });
        attr.write(path);
    }

    {
        Graph graph(n);
        auto attr = graph.nodeAttributes().attach<double>(name);
        attr.read(path);
        graph.forNodes([&](node u) { EXPECT_EQ(values[u], attr[u]); });
    }

    std::remove(path.c_str());
}

TEST_F(AttributeTest, testSetGetOnNonExistingNodes) {

    Graph graph(5);

    graph.removeNode(3);

    auto intAttr = graph.nodeAttributes().attach<int>("some int attribute");

    auto readAtIndex = [&intAttr](node n) -> int { return intAttr[n]; };

    EXPECT_FALSE(graph.hasNode(5));
    EXPECT_THROW(intAttr.set(5, 5), std::runtime_error);
    EXPECT_THROW(intAttr.get(5), std::runtime_error);
    EXPECT_THROW(intAttr[5] = 5, std::runtime_error);
    // trigger read access by int conversion
    EXPECT_THROW(readAtIndex(5), std::runtime_error);

    EXPECT_FALSE(graph.hasNode(3));
    EXPECT_THROW(intAttr.set(3, 3), std::runtime_error);
    EXPECT_THROW(intAttr.get(3), std::runtime_error);
    EXPECT_THROW(intAttr[3] = 3, std::runtime_error);
    // trigger read access by int conversion
    EXPECT_THROW(readAtIndex(3), std::runtime_error);
}

TEST_F(AttributeTest, testAttachDetachAttach) {

    Graph graph(15);

    auto intAttr = graph.nodeAttributes().attach<int>("some int attribute");

    intAttr[3] = 33;
    EXPECT_EQ(33, intAttr[3]);

    graph.nodeAttributes().detach("some int attribute");
    EXPECT_THROW(intAttr.get(3), std::runtime_error);

    intAttr = graph.nodeAttributes().attach<int>("some new int attribute");

    intAttr[3] = 33;
    EXPECT_EQ(intAttr[3], 33);
}

TEST_F(AttributeTest, testDoubleAttach) {

    Graph graph(15);

    auto intAttr = graph.nodeAttributes().attach<int>("some int attribute");
    auto dblAttr = graph.nodeAttributes().attach<double>("some double attribute");

    EXPECT_THROW(graph.nodeAttributes().attach<int>("some int attribute"), std::runtime_error);
}

} // namespace NetworKit
