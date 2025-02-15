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
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class AttributeGTest : public testing::Test {};

template <class Type>
class AttributeTypeGTest : public testing::Test {
protected:
    using AttrType = Type;
};

using AttrTypes = testing::Types<int, double, std::string>;
TYPED_TEST_SUITE(AttributeTypeGTest, AttrTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(AttributeTypeGTest, testReadWrite) {
    using AttrType = typename TestFixture::AttrType;

    const int n = 10;
    const std::string path = "output/attribute.txt";
    const std::string name = "attribute";

    auto attrValue = [](node u, node v) {
        if constexpr (std::is_same_v<AttrType, int>)
            return AttrType{static_cast<int>(u + v)};
        if constexpr (std::is_same_v<AttrType, double>)
            return AttrType{u + v + 0.2};
        if constexpr (std::is_same_v<AttrType, std::string>)
            return AttrType{"AttributeGTest with spaces and, and special symbols !: and \t tabs"};
    };

    {
        Graph graph(n, false, false, true);
        graph.addEdge(0, 1);
        graph.addEdge(1, 2);
        graph.addEdge(2, 3);
        auto attr = graph.edgeAttributes().attach<AttrType>(name);
        graph.forEdges([&](node u, node v) { attr.set2(u, v, attrValue(u, v)); });
        attr.write(path);
    }

    {
        Graph graph(n, false, false, true);
        graph.addEdge(0, 1);
        graph.addEdge(1, 2);
        graph.addEdge(2, 3);
        auto attr = graph.edgeAttributes().attach<AttrType>(name);
        attr.read(path);
        graph.forEdges([&](node u, node v) { EXPECT_EQ(attrValue(u, v), attr.get2(u, v)); });
    }

    std::remove(path.c_str());
}

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

TEST_F(AttributeGTest, testConstGetNodeAttribute) {
    const Graph graph = []() {
        Graph graph(10);

        auto intAttr = graph.nodeAttributes().attach<int>("some int attribute");
        intAttr.set(0, 1);

        return graph;
    }();

    auto constAttr = graph.nodeAttributes().get<int>("some int attribute");
    EXPECT_EQ(constAttr[0], 1);
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

TEST_F(AttributeGTest, testConstGetEdgeAttribute) {
    const Graph graph = []() {
        Graph graph(10);
        graph.addEdge(0, 1);
        graph.indexEdges();

        auto edgeAttr = graph.edgeAttributes().attach<double>("some edge attribute");
        edgeAttr(0, 1) = 3;

        return graph;
    }();

    auto constEdgeAttr = graph.edgeAttributes().get<double>("some edge attribute");
    EXPECT_EQ(constEdgeAttr(0, 1), 3);
}

/// COMBINED ATTRIBUTE TESTS ///

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

TEST_F(AttributeGTest, testAttributeCopy) {

    Graph graph1(5);
    graph1.addEdge(0, 1);
    graph1.indexEdges();

    auto edgeAttr1 = graph1.edgeAttributes().attach<int>("int edge attr");
    edgeAttr1(0, 1) = 3;

    // make copy of graph (and its attributes)
    Graph graph2{graph1};

    // modify attribute of graph1
    edgeAttr1(0, 1) = 4;

    auto edgeAttr2 = graph2.edgeAttributes().get<int>("int edge attr");

    EXPECT_EQ(edgeAttr2(0, 1), 3);

    // modify graph1
    graph1.removeEdge(0, 1);

    EXPECT_EQ(edgeAttr2(0, 1), 3);
}

TEST_F(AttributeGTest, testGraphModeChange) {
    Graph dgraph(5, false, true);
    dgraph.indexEdges();
    dgraph.addEdge(0, 1);
    auto edgeAttr = dgraph.edgeAttributes().attach<int>("int edge attr");
    edgeAttr(0, 1) = 3;
    auto nodeAttr = dgraph.nodeAttributes().attach<int>("int node attr");
    nodeAttr[2] = 4;

    auto wdgraph = GraphTools::toWeighted(dgraph);
    auto wNodeAttr = wdgraph.nodeAttributes().get<int>("int node attr");
    EXPECT_EQ(wNodeAttr[2], 4);
    auto wEdgeAttr = wdgraph.edgeAttributes().get<int>("int edge attr");
    EXPECT_EQ(wEdgeAttr(0, 1), 3);

    auto graph3 = GraphTools::toUnweighted(wdgraph);
    auto wNodeAttr3 = graph3.nodeAttributes().get<int>("int node attr");
    EXPECT_EQ(wNodeAttr3[2], 4);
    auto wEdgeAttr3 = graph3.edgeAttributes().get<int>("int edge attr");
    EXPECT_EQ(wEdgeAttr3(0, 1), 3);
}

TEST_F(AttributeGTest, testMoveGraph) {
    Graph graph1(3);
    graph1.indexEdges();
    graph1.addEdge(0, 1);
    auto attr1 = graph1.edgeAttributes().attach<int>("attr");
    attr1(0, 1) = 1;
    Graph graph2(std::move(graph1));
    // make sure that attr2 does not point to graph1 which was moved and may not be accessed anymore
    auto attr2 = graph2.edgeAttributes().get<int>("attr");
    EXPECT_EQ(attr2(0, 1), 1);
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

TEST_F(AttributeGTest, deleteOtherNode) {
    Graph graph(3, false, true, true);

    graph.addEdge(0, 1);
    graph.addEdge(1, 2);

    auto nodeAttr = graph.nodeAttributes().attach<node>("nodeAttr");
    nodeAttr[0] = 0;
    nodeAttr[1] = 1;
    nodeAttr[2] = 2;

    std::cout << nodeAttr[0] << "\n";
    std::cout << nodeAttr[1] << "\n";
    std::cout << nodeAttr[2] << "\n";

    graph.removeNode(1);

    EXPECT_EQ(nodeAttr[0], 0);
    EXPECT_THROW(nodeAttr[1], std::runtime_error);
    EXPECT_EQ(nodeAttr[2], 2);
    std::cout << nodeAttr[2] << "\n";

    graph.removeNode(0);

    EXPECT_THROW(nodeAttr[0], std::runtime_error);
    EXPECT_THROW(nodeAttr[1], std::runtime_error);
    EXPECT_EQ(nodeAttr[2], 2);
}

TEST_F(AttributeGTest, testDeleteEdges) {
    NetworKit::Graph graph(3);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);

    graph.indexEdges();
    auto attr = graph.edgeAttributes().attach<double>("attr");

    attr(0, 1) = 1.0;
    attr(0, 2) = 2.0;

    graph.removeEdge(2, 0);

    EXPECT_EQ(attr(0, 1), 1.0);
}

TEST_F(AttributeGTest, assignAttributeToAttribute) {
    Graph graph(3);
    auto attr1 = graph.nodeAttributes().attach<int>("attr1");
    auto attr2 = graph.nodeAttributes().attach<int>("attr2");

    attr1[0] = 1;
    attr2[0] = 3;
    attr2[0] = attr1[0]; // operator=(IndexProxy&&)
    EXPECT_EQ(attr2[0], 1);

    attr1[1] = 4;
    auto v1 = attr1[1];
    EXPECT_EQ(v1, 4);
    auto v2 = attr2[1];
    v2 = v1; // operator=(const IndexProxy&)
    EXPECT_EQ(v2, 4);
    EXPECT_EQ(attr2[1], 4);
}

TEST_F(AttributeGTest, testMaintainCompactEdges) {
    Graph graph(3);
    graph.indexEdges();
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);

    graph.setMaintainCompactEdges();

    auto attr = graph.edgeAttributes().attach<int>("attr");
    attr(0, 1) = 0;
    attr(0, 2) = 1;

    graph.removeEdge(0, 1);

    EXPECT_EQ(attr(0, 2), 1);
}

TEST_F(AttributeGTest, testMaintainSortedEdges) {
    Graph graph(3);
    graph.indexEdges();
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);

    graph.setKeepEdgesSorted();

    auto attr = graph.edgeAttributes().attach<int>("attr");
    attr(0, 1) = 0;
    attr(0, 2) = 1;

    graph.removeEdge(0, 1);

    EXPECT_EQ(attr(0, 2), 1);
}

TEST_F(AttributeGTest, testMaintainSortedAndCompactEdges) {
    Graph graph(3);
    graph.indexEdges();
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);

    graph.setMaintainCompactEdges();
    graph.setKeepEdgesSorted();

    auto attr = graph.edgeAttributes().attach<int>("attr");
    attr(0, 1) = 0;
    attr(0, 2) = 1;

    graph.removeEdge(0, 1);

    EXPECT_EQ(attr(0, 2), 1);
}

} // namespace NetworKit
