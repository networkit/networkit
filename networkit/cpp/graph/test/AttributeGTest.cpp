/*
 * AttributeGTest.cpp
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

class AttributeGTest : public testing::Test {};

/// NODE ATTRIBUTE TESTS ///

TEST_F(AttributeGTest, testNodeAttributeSetGetOnExistingNodes) {

    Graph graph(15);

    const std::string name = "some int attribute";
    auto intAttr = graph.nodeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.nodeAttributes().find(name)->second->getName(), name);

    // set / get
    graph.forNodes([&](node n) { intAttr.set(n, int(n)); });
    EXPECT_EQ(intAttr.size(), graph.numberOfNodes());

    graph.forNodes([&](node n) { EXPECT_EQ(intAttr.get(n), int(n)); });
}

TEST_F(AttributeGTest, testNodeAttributeSetGetOnExistingNodesByIndexProxy) {

    Graph graph(15);

    const std::string name = "some int attribute";
    auto intAttr = graph.nodeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.nodeAttributes().find(name)->second->getName(), name);

    // by index proxy
    graph.forNodes([&](node n) { intAttr[n] = int(n); });
    EXPECT_EQ(intAttr.size(), graph.numberOfNodes());

    graph.forNodes([&](node n) { EXPECT_EQ(intAttr[n], int(n)); });
}

TEST_F(AttributeGTest, testNodeAttributeIteratorOnExistingNodes) {

    Graph graph(15);

    const std::string name = "some int attribute";
    auto intAttr = graph.nodeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.nodeAttributes().find(name)->second->getName(), name);

    // by index proxy
    graph.forNodes([&](node n) { intAttr[n] = int(n); });

    // Test iterator
    node u = 0;
    int att = 0;
    for (auto pair : intAttr) {
        EXPECT_EQ(pair.first, u);
        EXPECT_EQ(pair.second, att);
        ++u, ++att;
    }
}

TEST_F(AttributeGTest, testNodeAttributeReadWrite) {
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

/// EDGE ATTRIBUTE TESTS ///

TEST_F(AttributeGTest, testEdgeAttributeSetGetOnExistingEdges) {

    Graph graph(15, false, false, true);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);

    const std::string name = "some int attribute";
    auto intAttr = graph.edgeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.edgeAttributes().find(name)->second->getName(), name);

    // set / get
    graph.forEdges([&](node u, node v) { intAttr.set2(u, v, int(u) + int(v)); });
    EXPECT_EQ(intAttr.size(), graph.numberOfEdges());

    graph.forEdges([&](node u, node v) { EXPECT_EQ(intAttr.get2(u, v), int(u) + int(v)); });
}

TEST_F(AttributeGTest, testEdgeAttributeSetGetOnExistingEdgesByIndexProxy) {

    Graph graph(15, false, false, true);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);

    const std::string name = "some int attribute";
    auto intAttr = graph.edgeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.edgeAttributes().find(name)->second->getName(), name);

    // by index proxy
    graph.forEdges([&](node u, node v) { intAttr(u, v) = int(u) + int(v); });
    EXPECT_EQ(intAttr.size(), graph.numberOfEdges());

    graph.forEdges([&](node u, node v) { EXPECT_EQ(intAttr(u, v), int(u) + int(v)); });
}

TEST_F(AttributeGTest, testEdgeAttributeSetGetOnExistingEdgesByIndexedEdge) {

    Graph graph(15, false, false, true);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);

    const std::string name = "some int attribute";
    auto intAttr = graph.edgeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.edgeAttributes().find(name)->second->getName(), name);

    // by index proxy
    graph.forEdges([&](node u, node v) {
        auto edgeIndex = graph.edgeId(u, v);
        intAttr[edgeIndex] = int(u) + int(v);
    });
    EXPECT_EQ(intAttr.size(), graph.numberOfEdges());

    graph.forEdges([&](node u, node v) {
        auto edgeIndex = graph.edgeId(u, v);
        EXPECT_EQ(intAttr[edgeIndex], int(u) + int(v));
    });
}

TEST_F(AttributeGTest, testEdgeAttributeIteratorOnExistingEdges) {

    Graph graph(15, false, false, true);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);

    const std::string name = "some int attribute";
    auto intAttr = graph.edgeAttributes().attach<int>(name);
    EXPECT_EQ(intAttr.size(), 0u);
    EXPECT_EQ(graph.edgeAttributes().find(name)->second->getName(), name);

    // by index proxy
    graph.forEdges([&](node u, node v) { intAttr(u, v) = u + v; });

    // Test iterator
    node u = 0;
    int att = 0;
    for (auto pair : intAttr) {
        EXPECT_EQ(pair.first, u);
        auto edge = graph.edgeById(u);
        EXPECT_EQ(pair.second, edge.first + edge.second);
        ++u, ++att;
    }
}

TEST_F(AttributeGTest, testEdgeAttributeReadWrite) {
    const int n = 10;
    const std::string path = "output/attribute.txt";
    const std::string name = "attribute";

    {
        Graph graph(n, false, false, true);
        graph.addEdge(0, 1);
        graph.addEdge(1, 2);
        graph.addEdge(2, 3);
        auto attr = graph.edgeAttributes().attach<double>(name);
        graph.forEdges([&](node u, node v) { attr.set2(u, v, u + v); });
        attr.write(path);
    }

    {
        Graph graph(n, false, false, true);
        graph.addEdge(0, 1);
        graph.addEdge(1, 2);
        graph.addEdge(2, 3);
        auto attr = graph.edgeAttributes().attach<double>(name);
        attr.read(path);
        graph.forEdges([&](node u, node v) { EXPECT_EQ(u + v, attr.get2(u, v)); });
    }

    std::remove(path.c_str());
}

/// COMBINED ATTRIBUTE TESTS ///

TEST_F(AttributeGTest, testAttributeSetGetOnNonExistingNodes) {

    auto tester = [&](auto &attributeMap) {
        auto intAttr = attributeMap.template attach<int>("some int attribute");

        auto readAtIndex = [&intAttr](node n) -> int { return intAttr[n]; };

        EXPECT_THROW(intAttr.set(5, 5), std::runtime_error);
        EXPECT_THROW(intAttr.get(5), std::runtime_error);
        EXPECT_THROW(intAttr[5] = 5, std::runtime_error);
        // trigger read access by int conversion
        EXPECT_THROW(readAtIndex(5), std::runtime_error);

        EXPECT_THROW(intAttr.set(3, 3), std::runtime_error);
        EXPECT_THROW(intAttr.get(3), std::runtime_error);
        EXPECT_THROW(intAttr[3] = 3, std::runtime_error);
        // trigger read access by int conversion
        EXPECT_THROW(readAtIndex(3), std::runtime_error);

        attributeMap.detach("some int attribute");
    };

    Graph graph(5, false, false, true);
    graph.removeNode(3);

    tester(graph.nodeAttributes());
    tester(graph.edgeAttributes());
}

TEST_F(AttributeGTest, testAttributeAttachDetachAttach) {

    auto tester = [&](auto &attributeMap) {
        auto intAttr = attributeMap.template attach<int>("some int attribute");
        intAttr[3] = 33;
        EXPECT_EQ(33, intAttr[3]);

        attributeMap.detach("some int attribute");
        EXPECT_THROW(intAttr.get(3), std::runtime_error);

        intAttr = attributeMap.template attach<int>("some new int attribute");

        intAttr[3] = 33;
        EXPECT_EQ(intAttr[3], 33);

        attributeMap.detach("some new int attribute");
    };

    Graph graph(15, false, false, true);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);

    tester(graph.nodeAttributes());
    tester(graph.edgeAttributes());
}

TEST_F(AttributeGTest, testAttributeDoubleAttach) {

    auto tester = [&](auto &attributeMap) {
        auto intAttr = attributeMap.template attach<int>("some int attribute");
        auto dblAttr = attributeMap.template attach<double>("some double attribute");

        EXPECT_THROW(attributeMap.template attach<int>("some int attribute"), std::runtime_error);

        attributeMap.detach("some int attribute");
        attributeMap.detach("some double attribute");
    };

    Graph graph(15, false, false, true);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);

    tester(graph.nodeAttributes());
    tester(graph.edgeAttributes());
}

TEST_F(AttributeGTest, testInvalidate) {

    Graph graph(10);

    auto intAttr = graph.nodeAttributes().attach<int>("some int attribute");
    intAttr.set(0, 1);
    EXPECT_EQ(intAttr.size(), 1);
    // deleting the node with the attached attribute
    graph.removeNode(0);
    // the invalidate method should have reset the intAttr
    EXPECT_EQ(intAttr.size(), 0);
}

TEST_F(AttributeGTest, testDefaultGet) {

    Graph graph(10);

    auto intAttr = graph.nodeAttributes().attach<int>("some int attribute");
    intAttr.set(0, 1);
    // trying to get attribute at invalid index (2), leads to default return
    EXPECT_EQ(intAttr.get(2, 0), 0);
}

} // namespace NetworKit
