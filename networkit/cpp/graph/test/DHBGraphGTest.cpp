#include <gtest/gtest.h>
#include <networkit/graph/DHBGraph.hpp>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>

#include <networkit/io/METISGraphReader.hpp>

using namespace NetworKit;

class DHBGraphGTest : public ::testing::TestWithParam<std::tuple<bool, bool>> {
public:
    virtual void SetUp();

protected:
    DHBGraph Ghouse;
    std::vector<std::pair<node, node>> houseEdgesOut;
    std::vector<std::vector<edgeweight>> Ahouse;
    count n_house;
    count m_house;

    bool isGraph() const { return !isWeighted() && !isDirected(); }
    bool isWeightedGraph() const { return isWeighted() && !isDirected(); }
    bool isDirectedGraph() const { return !isWeighted() && isDirected(); }
    bool isWeightedDirectedGraph() const { return isWeighted() && isDirected(); }

    bool isWeighted() const;
    bool isDirected() const;
    DHBGraph createGraph(count n = 0) const;
    DHBGraph createGraph(count n, count m) const;
    count countSelfLoopsManually(const DHBGraph &G);
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, DHBGraphGTest,
                         testing::Values(std::make_tuple(false, false),
                                         std::make_tuple(true, false), std::make_tuple(false, true),
                                         std::make_tuple(true, true)));

bool DHBGraphGTest::isWeighted() const {
    return std::get<0>(GetParam());
}
bool DHBGraphGTest::isDirected() const {
    return std::get<1>(GetParam());
}

DHBGraph DHBGraphGTest::createGraph(count n) const {
    bool weighted, directed;
    std::tie(weighted, directed) = GetParam();
    DHBGraph G(n, weighted, directed);
    return G;
}

DHBGraph DHBGraphGTest::createGraph(count n, count m) const {
    auto G = createGraph(n);
    while (G.numberOfEdges() < m) {
        const auto u = Aux::Random::index(n);
        const auto v = Aux::Random::index(n);
        if (u == v)
            continue;
        if (G.hasEdge(u, v))
            continue;

        const auto p = Aux::Random::probability();
        G.addEdge(u, v, p);
    }
    return G;
}

count DHBGraphGTest::countSelfLoopsManually(const DHBGraph &G) {
    count c = 0;
    G.forEdges([&](node u, node v) {
        if (u == v) {
            c += 1;
        }
    });
    return c;
}

void DHBGraphGTest::SetUp() {
    /*
     *    0
     *   . \
     *  /   \
     * /     .
     * 1 <-- 2
     * ^ \  .|
     * |  \/ |
     * | / \ |
     * |/   ..
     * 3 <-- 4
     *
     * move you pen from node to node:
     * 3 -> 1 -> 0 -> 2 -> 1 -> 4 -> 3 -> 2 -> 4
     */
    n_house = 5;
    m_house = 8;

    Ghouse = createGraph(5);
    houseEdgesOut = {{3, 1}, {1, 0}, {0, 2}, {2, 1}, {1, 4}, {4, 3}, {3, 2}, {2, 4}};
    Ahouse = {n_house, std::vector<edgeweight>(n_house, 0.0)};
    edgeweight ew = 1.0;
    for (auto &e : houseEdgesOut) {
        node u = e.first;
        node v = e.second;
        Ghouse.addEdge(u, v, ew);

        Ahouse[u][v] = ew;

        if (!Ghouse.isDirected()) {
            Ahouse[v][u] = ew;
        }

        if (Ghouse.isWeighted()) {
            ew += 1.0;
        }
    }
}

/** CONSTRUCTORS **/

TEST(DHBGraphGTest, testDefConstructorWithUndirIndex) {
    // Test indexed + undirected graph
    DHBGraph GUndir(3, false, false, true);
    GUndir.addEdge(0, 1);
    EXPECT_EQ(GUndir.edgeId(0, 1), 0);
    EXPECT_EQ(GUndir.edgeId(1, 0), 0);
}

TEST(DHBGraphGTest, testDefConstructorWithDirIndex) {
    // Test indexed + directed graph
    DHBGraph GDir(3, false, true, true);
    GDir.addEdge(0, 1);
    EXPECT_EQ(GDir.edgeId(0, 1), 0);
}

TEST_P(DHBGraphGTest, testCopyConstructorWithIndexedEdgeIds) {
    DHBGraph G(3, false, false, true);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges(true);

    DHBGraph GCopy(G, isWeighted(), isDirected(), true);
    EXPECT_TRUE(GCopy.hasEdgeIds());
    EXPECT_TRUE(GCopy.hasEdge(0, 1));
    EXPECT_TRUE(GCopy.hasEdge(1, 2));
}

TEST_P(DHBGraphGTest, testConstructorWithInitializerListWeightedEdges) {
    // Create a graph using an initializer list of WeightedEdges
    DHBGraph G({{0, 1, 1.5}, {1, 2, 2.5}, {2, 3, 3.5}, {3, 0, 4.5}});

    // Verify the graph properties
    ASSERT_EQ(4u, G.numberOfNodes());
    ASSERT_EQ(4u, G.numberOfEdges());

    // Check if the edges exist with the correct weights
    ASSERT_TRUE(G.hasEdge(0, 1));
    ASSERT_TRUE(G.hasEdge(1, 2));
    ASSERT_TRUE(G.hasEdge(2, 3));
    ASSERT_TRUE(G.hasEdge(3, 0));

    if (G.isWeighted()) {
        ASSERT_EQ(1.5, G.weight(0, 1));
        ASSERT_EQ(2.5, G.weight(1, 2));
        ASSERT_EQ(3.5, G.weight(2, 3));
        ASSERT_EQ(4.5, G.weight(3, 0));
    } else {
        ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
        ASSERT_EQ(defaultEdgeWeight, G.weight(1, 2));
        ASSERT_EQ(defaultEdgeWeight, G.weight(2, 3));
        ASSERT_EQ(defaultEdgeWeight, G.weight(3, 0));
    }

    if (!G.isDirected()) {
        // Check symmetry for undirected graphs
        ASSERT_TRUE(G.hasEdge(1, 0));
        ASSERT_TRUE(G.hasEdge(2, 1));
        ASSERT_TRUE(G.hasEdge(3, 2));
        ASSERT_TRUE(G.hasEdge(0, 3));

        if (G.isWeighted()) {
            ASSERT_EQ(1.5, G.weight(1, 0));
            ASSERT_EQ(2.5, G.weight(2, 1));
            ASSERT_EQ(3.5, G.weight(3, 2));
            ASSERT_EQ(4.5, G.weight(0, 3));
        } else {
            ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
            ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
            ASSERT_EQ(defaultEdgeWeight, G.weight(3, 2));
            ASSERT_EQ(defaultEdgeWeight, G.weight(0, 3));
        }
    }
}

TEST_P(DHBGraphGTest, testCopyConstructor) {
    // G, GW, D, DW - the copy
    // This - the origin
    DHBGraph G = DHBGraph(this->Ghouse, false, false);
    DHBGraph GW = DHBGraph(this->Ghouse, true, false);
    DHBGraph D = DHBGraph(this->Ghouse, false, true);
    DHBGraph DW = DHBGraph(this->Ghouse, true, true);

    ASSERT_FALSE(G.isWeighted());
    ASSERT_FALSE(G.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), G.numberOfNodes());

    // copying from undirected graph to directed graph
    ASSERT_EQ(this->Ghouse.numberOfEdges(), G.numberOfEdges());

    ASSERT_TRUE(GW.isWeighted());
    ASSERT_FALSE(GW.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), GW.numberOfNodes());

    ASSERT_FALSE(D.isWeighted());
    ASSERT_TRUE(D.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), D.numberOfNodes());

    // if orig is not directed, but copy to a directed graph, edge number should double.
    ASSERT_EQ(this->Ghouse.numberOfEdges(), D.numberOfEdges());

    ASSERT_TRUE(DW.isWeighted());
    ASSERT_TRUE(DW.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), DW.numberOfNodes());
    ASSERT_EQ(this->Ghouse.numberOfEdges(), DW.numberOfEdges());
    this->Ghouse.forNodes([&](node v) {
        count d = this->Ghouse.degree(v);
        // count dUndirected = isDirected() ? d + this->Ghouse.degreeIn(v) : d;

        // Comment out following evaluation, since when copying from an directed graph to an
        // undirected graph, the degree of a vertex may change.

        // ASSERT_EQ(dUndirected, G.degree(v));
        // ASSERT_EQ(dUndirected, GW.degree(v));
        ASSERT_EQ(d, D.degree(v));
        ASSERT_EQ(d, DW.degree(v));
    });

    // if Ghouse was directed we should have an exact copy of it, but if it was
    // undirected we should have edges in both directions
    count m = 0;
    G.forEdges([&](node u, node v) {
        ASSERT_TRUE(G.hasEdge(v, u));
        ASSERT_EQ(defaultEdgeWeight, G.weight(v, u));
        ASSERT_EQ(defaultEdgeWeight, G.weight(u, v));

        auto e = std::make_pair(u, v);
        bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                     != this->houseEdgesOut.end();
        if (!found) {
            e = std::make_pair(v, u);
            found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                    != this->houseEdgesOut.end();
        }
        ASSERT_TRUE(found);
        m++;
    });
    ASSERT_EQ(8u, m);

    m = 0;
    GW.forEdges([&](node u, node v) {
        ASSERT_TRUE(GW.hasEdge(v, u));
        ASSERT_EQ(GW.weight(u, v), GW.weight(v, u));

        auto e = std::make_pair(u, v);
        bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                     != this->houseEdgesOut.end();
        if (!found) {
            e = std::make_pair(v, u);
            found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                    != this->houseEdgesOut.end();
            ASSERT_EQ(this->Ghouse.weight(v, u), GW.weight(v, u));
        } else {
            ASSERT_EQ(this->Ghouse.weight(u, v), GW.weight(u, v));
        }
        ASSERT_TRUE(found);
        m++;
    });
    ASSERT_EQ(8u, m);

    m = 0;
    D.forEdges([&](node u, node v) {
        ASSERT_EQ(defaultEdgeWeight, D.weight(u, v));

        auto e = std::make_pair(u, v);
        bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                     != this->houseEdgesOut.end();
        if (!this->Ghouse.isDirected()) {
            ASSERT_TRUE(D.hasEdge(v, u));

            e = std::make_pair(v, u);
            found = found
                    || (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                        != this->houseEdgesOut.end());
        } else {
            ASSERT_FALSE(D.hasEdge(v, u));
        }
        ASSERT_TRUE(found);
        m++;
    });
    count m_expected = isDirected() ? 8 : 16;
    ASSERT_EQ(m_expected, m);

    m = 0;
    DW.forEdges([&](node u, node v) {
        ASSERT_EQ(this->Ghouse.weight(u, v), DW.weight(u, v));

        auto e = std::make_pair(u, v);
        bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                     != this->houseEdgesOut.end();
        if (!this->Ghouse.isDirected()) {
            ASSERT_TRUE(DW.hasEdge(v, u));
            e = std::make_pair(v, u);
            found = found
                    || (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e)
                        != this->houseEdgesOut.end());
        } else {
            ASSERT_FALSE(DW.hasEdge(v, u));
        }
        ASSERT_TRUE(found);
        m++;
    });
    m_expected = isDirected() ? 8 : 16;
    ASSERT_EQ(m_expected, m);
}

/** NODE MODIFIERS **/

TEST_P(DHBGraphGTest, testAddNode) {
    DHBGraph G = createGraph();

    ASSERT_FALSE(G.hasNode(0));
    ASSERT_FALSE(G.hasNode(1));
    ASSERT_EQ(0u, G.numberOfNodes());

    node node_0 = G.addNode();
    ASSERT_EQ(0, node_0);
    ASSERT_TRUE(G.hasNode(0));
    ASSERT_FALSE(G.hasNode(1));
    ASSERT_EQ(1u, G.numberOfNodes());

    DHBGraph G2 = createGraph(2);
    ASSERT_TRUE(G2.hasNode(0));
    ASSERT_TRUE(G2.hasNode(1));
    ASSERT_FALSE(G2.hasNode(2));
    ASSERT_EQ(2u, G2.numberOfNodes());

    node node_2 = G2.addNode();
    node node_3 = G2.addNode();
    ASSERT_EQ(2, node_2);
    ASSERT_EQ(3, node_3);
    ASSERT_TRUE(G2.hasNode(2));
    ASSERT_TRUE(G2.hasNode(3));
    ASSERT_FALSE(G2.hasNode(4));
    ASSERT_EQ(4u, G2.numberOfNodes());
}

TEST_P(DHBGraphGTest, testAddNodes) {
    auto G = createGraph(5);
    node add_5 = G.addNodes(5);

    ASSERT_EQ(9, add_5);

    ASSERT_EQ(G.numberOfNodes(), 10);
    ASSERT_TRUE(G.hasNode(9));

    G.addNodes(90);

    ASSERT_EQ(G.numberOfNodes(), 100);
    ASSERT_TRUE(G.hasNode(99));
}

TEST_P(DHBGraphGTest, testHasNode) {
    DHBGraph G = createGraph(5);

    ASSERT_TRUE(G.hasNode(0));
    ASSERT_TRUE(G.hasNode(1));
    ASSERT_TRUE(G.hasNode(2));
    ASSERT_TRUE(G.hasNode(3));
    ASSERT_TRUE(G.hasNode(4));
    ASSERT_FALSE(G.hasNode(5));
    ASSERT_FALSE(G.hasNode(6));
}

TEST_P(DHBGraphGTest, testRemoveNodeDHB) {
    this->Ghouse.removeNode(0);
    ASSERT_TRUE(this->Ghouse.isIsolated(0));
    this->Ghouse.removeNode(1);
    ASSERT_TRUE(this->Ghouse.isIsolated(1));
    this->Ghouse.removeNode(2);
    ASSERT_TRUE(this->Ghouse.isIsolated(2));
    this->Ghouse.removeNode(3);
    ASSERT_TRUE(this->Ghouse.isIsolated(3));
    this->Ghouse.removeNode(4);
    ASSERT_TRUE(this->Ghouse.isIsolated(4));
}

/** NODE PROPERTIES **/

TEST_P(DHBGraphGTest, testDegree) {
    if (isDirected()) {
        ASSERT_EQ(1u, this->Ghouse.degree(0));
        ASSERT_EQ(2u, this->Ghouse.degree(1));
        ASSERT_EQ(2u, this->Ghouse.degree(2));
        ASSERT_EQ(2u, this->Ghouse.degree(3));
        ASSERT_EQ(1u, this->Ghouse.degree(4));
    } else {
        ASSERT_EQ(2u, this->Ghouse.degree(0));
        ASSERT_EQ(4u, this->Ghouse.degree(1));
        ASSERT_EQ(4u, this->Ghouse.degree(2));
        ASSERT_EQ(3u, this->Ghouse.degree(3));
        ASSERT_EQ(3u, this->Ghouse.degree(4));
    }
}

TEST_P(DHBGraphGTest, testDegreeIn) {
    if (isDirected()) {
        ASSERT_EQ(1u, this->Ghouse.degreeIn(0));
        ASSERT_EQ(2u, this->Ghouse.degreeIn(1));
        ASSERT_EQ(2u, this->Ghouse.degreeIn(2));
        ASSERT_EQ(1u, this->Ghouse.degreeIn(3));
        ASSERT_EQ(2u, this->Ghouse.degreeIn(4));
    } else {
        ASSERT_EQ(2u, this->Ghouse.degreeIn(0));
        ASSERT_EQ(4u, this->Ghouse.degreeIn(1));
        ASSERT_EQ(4u, this->Ghouse.degreeIn(2));
        ASSERT_EQ(3u, this->Ghouse.degreeIn(3));
        ASSERT_EQ(3u, this->Ghouse.degreeIn(4));
    }
}

TEST_P(DHBGraphGTest, testDegreeOut) {
    if (isDirected()) {
        ASSERT_EQ(1u, this->Ghouse.degreeOut(0));
        ASSERT_EQ(2u, this->Ghouse.degreeOut(1));
        ASSERT_EQ(2u, this->Ghouse.degreeOut(2));
        ASSERT_EQ(2u, this->Ghouse.degreeOut(3));
        ASSERT_EQ(1u, this->Ghouse.degreeOut(4));
    } else {
        ASSERT_EQ(2u, this->Ghouse.degreeOut(0));
        ASSERT_EQ(4u, this->Ghouse.degreeOut(1));
        ASSERT_EQ(4u, this->Ghouse.degreeOut(2));
        ASSERT_EQ(3u, this->Ghouse.degreeOut(3));
        ASSERT_EQ(3u, this->Ghouse.degreeOut(4));
    }
}

TEST_P(DHBGraphGTest, testIsIsolated) {
    ASSERT_FALSE(this->Ghouse.isIsolated(0));
    ASSERT_FALSE(this->Ghouse.isIsolated(1));
    ASSERT_FALSE(this->Ghouse.isIsolated(2));
    ASSERT_FALSE(this->Ghouse.isIsolated(3));
    ASSERT_FALSE(this->Ghouse.isIsolated(4));

    this->Ghouse.addNode();
    ASSERT_TRUE(this->Ghouse.isIsolated(5));

    this->Ghouse.removeEdge(1, 0);
    ASSERT_FALSE(this->Ghouse.isIsolated(0));

    this->Ghouse.removeEdge(0, 2);
    ASSERT_TRUE(this->Ghouse.isIsolated(0));

    this->Ghouse.addEdge(1, 0);
    ASSERT_FALSE(this->Ghouse.isIsolated(0));
}

TEST_P(DHBGraphGTest, testWeightedDegree) {
    // add self-loop
    this->Ghouse.addEdge(2, 2, 0.75);

    if (isGraph()) {
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(5 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(5.0, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(12.0, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(22.75, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(14.0, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(19.0, this->Ghouse.weightedDegree(4));
    }

    if (isDirectedGraph()) {
        // only count outgoing edges
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
    }

    if (isWeightedDirectedGraph()) {
        // only sum weight of outgoing edges
        ASSERT_EQ(3.0, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(7.0, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(12.75, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(8.0, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(6.0, this->Ghouse.weightedDegree(4));
    }
}

TEST_P(DHBGraphGTest, testWeightedDegree2) {
    // add self-loop
    this->Ghouse.addEdge(2, 2, 0.75);

    if (isGraph()) {
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(6 * defaultEdgeWeight, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(4, true));
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(5.0, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(12.0, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(23.5, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(14.0, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(19.0, this->Ghouse.weightedDegree(4, true));
    }

    if (isDirectedGraph()) {
        // only count outgoing edges
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(4, true));
    }

    if (isWeightedDirectedGraph()) {
        // only sum weight of outgoing edges
        ASSERT_EQ(3.0, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(7.0, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(13.5, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(8.0, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(6.0, this->Ghouse.weightedDegree(4, true));
    }
}

/** EDGE MODIFIERS **/

TEST_P(DHBGraphGTest, testIndexEdgesForce) {
    DHBGraph G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges(true);

    EXPECT_EQ(G.edgeId(0, 1), 0);
    EXPECT_EQ(G.edgeId(1, 2), 1);
    if (isGraph()) {
        EXPECT_EQ(G.edgeId(1, 0), 0);
        EXPECT_EQ(G.edgeId(2, 1), 1);
    }
}

TEST_P(DHBGraphGTest, testAddEdge) {
    DHBGraph G = createGraph(3);

    // DHBGraph without edges
    ASSERT_EQ(0u, G.numberOfEdges());
    ASSERT_FALSE(G.hasEdge(0, 2));
    ASSERT_FALSE(G.hasEdge(0, 1));
    ASSERT_FALSE(G.hasEdge(1, 2));
    ASSERT_FALSE(G.hasEdge(2, 2));
    ASSERT_EQ(nullWeight, G.weight(0, 2));
    ASSERT_EQ(nullWeight, G.weight(0, 1));
    ASSERT_EQ(nullWeight, G.weight(1, 2));
    ASSERT_EQ(nullWeight, G.weight(2, 2));

    // DHBGraph with 2 normal edges
    G.addEdge(0, 1, 4.51);
    G.addEdge(1, 2, 2.39);
    ASSERT_EQ(2u, G.numberOfEdges());
    ASSERT_FALSE(G.hasEdge(0, 2)); // was never added
    ASSERT_TRUE(G.hasEdge(0, 1));
    ASSERT_TRUE(G.hasEdge(1, 2));
    ASSERT_FALSE(G.hasEdge(2, 2)); // will be added later

    // check weights
    if (G.isWeighted()) {
        ASSERT_EQ(4.51, G.weight(0, 1));
        ASSERT_EQ(2.39, G.weight(1, 2));
    } else {
        ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
        ASSERT_EQ(defaultEdgeWeight, G.weight(1, 2));
    }

    if (G.isDirected()) {
        ASSERT_FALSE(G.hasEdge(1, 0));
        ASSERT_FALSE(G.hasEdge(2, 1));

        // add edge in the other direction
        // note: bidirectional edges are not supported, so both edges have different
        // weights
        G.addEdge(2, 1, 6.23);
        ASSERT_TRUE(G.hasEdge(2, 1));
        if (G.isWeighted()) {
            ASSERT_EQ(2.39, G.weight(1, 2));
            ASSERT_EQ(6.23, G.weight(2, 1));
        } else {
            ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
        }
    } else {
        ASSERT_TRUE(G.hasEdge(1, 0));
        ASSERT_TRUE(G.hasEdge(2, 1));
        if (G.isWeighted()) {
            ASSERT_EQ(4.51, G.weight(1, 0));
            ASSERT_EQ(2.39, G.weight(2, 1));
        } else {
            ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
            ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
        }
    }

    // add self loop
    G.addEdge(2, 2, 0.72);
    ASSERT_TRUE(G.hasEdge(2, 2));
    if (G.isWeighted()) {
        ASSERT_EQ(0.72, G.weight(2, 2));
    } else {
        ASSERT_EQ(defaultEdgeWeight, G.weight(2, 2));
    }
}

TEST_P(DHBGraphGTest, testAddEdges_weighted_edge_no_update) {
    DHBGraph G = createGraph(5);
    WeightedEdge e1(0, 2, defaultEdgeWeight);
    WeightedEdge e2(1, 2, defaultEdgeWeight);
    std::vector<WeightedEdge> edges;
    edges.push_back(e1);
    edges.push_back(e2);

    G.addEdges(std::move(edges), false);

    ASSERT_TRUE(G.hasEdge(0, 2));
    ASSERT_TRUE(G.hasEdge(1, 2));
}

TEST_P(DHBGraphGTest, testAddEdges_large_graph) {
    uint32_t constexpr num_vertices = 100'000;
    DHBGraph G = createGraph(num_vertices);

    float constexpr erdos_renyi_p = 0.0001;
    ErdosRenyiGenerator gen(num_vertices, erdos_renyi_p);
    Graph er_graph = gen.generate();

    std::vector<Edge> edges;
    er_graph.forEdges([&edges](NetworKit::node u, NetworKit::node v) {
        edges.push_back(Edge{u, v});
    });

    G.addEdges(std::move(edges), false);

    for (auto const &e : edges) {
        NetworKit::node const u = e.u;
        NetworKit::node const v = e.v;
        ASSERT_TRUE(G.hasEdge(u, v)) << "Edge (" << u << ", " << v << ") should exist.";
    }
}

TEST_P(DHBGraphGTest, testAddEdges_weighted_edge_insert_twice_no_update) {
    DHBGraph G = createGraph(5);
    WeightedEdge e1(0, 2, defaultEdgeWeight);
    WeightedEdge e2(1, 2, defaultEdgeWeight);
    G.addEdge(0, 2, 5.0);
    G.addEdge(1, 2, 5.0);

    NetworKit::edgeweight new_weight = 2.0f;
    e1.weight = new_weight;
    e2.weight = new_weight;

    std::vector<WeightedEdge> edges;
    edges.push_back(e1);
    edges.push_back(e2);

    G.addEdges(std::move(edges), false);

    ASSERT_TRUE(G.hasEdge(0, 2));
    ASSERT_TRUE(G.hasEdge(1, 2));

    ASSERT_FALSE(new_weight == G.weight(0, 2));
    ASSERT_FALSE(new_weight == G.weight(1, 2));
}

TEST_P(DHBGraphGTest, testAddEdges_weighted_edge_no_update_failed) {
    DHBGraph G = createGraph(5);
    WeightedEdge e1(0, 2, defaultEdgeWeight);
    WeightedEdge e2(1, 2, defaultEdgeWeight);
    G.addEdge(0, 2, 5.0);
    G.addEdge(1, 2, 5.0);
    std::vector<WeightedEdge> edges;
    edges.push_back(e1);
    edges.push_back(e2);
    G.addEdges(std::move(edges), false);
    ASSERT_TRUE(G.hasEdge(0, 2));
    ASSERT_TRUE(G.hasEdge(1, 2));
}

TEST_P(DHBGraphGTest, testAddEdges_edge_no_update) {
    DHBGraph G = createGraph(5);
    Edge e1(0, 2);
    Edge e2(1, 2);
    std::vector<Edge> edges;
    edges.push_back(e1);
    edges.push_back(e2);
    G.addEdges(std::move(edges), false);
    ASSERT_TRUE(G.hasEdge(0, 2));
    ASSERT_TRUE(G.hasEdge(1, 2));
}

TEST_P(DHBGraphGTest, testRemoveEdge) {
    double epsilon = 1e-6;
    DHBGraph G = createGraph(3);

    edgeweight ewBefore = G.totalEdgeWeight();

    G.addEdge(0, 1, 3.14);

    if (G.isWeighted()) {
        ASSERT_NEAR(ewBefore + 3.14, G.totalEdgeWeight(), epsilon);
    } else {
        ASSERT_NEAR(ewBefore + defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    G.addEdge(0, 0);

    ASSERT_EQ(2u, G.numberOfEdges());
    ASSERT_TRUE(G.hasEdge(0, 0));
    ASSERT_TRUE(G.hasEdge(0, 1));
    ASSERT_FALSE(G.hasEdge(2, 1));

    // test remove regular edge
    ewBefore = G.totalEdgeWeight();
    G.removeEdge(0, 1);
    if (G.isWeighted()) {
        ASSERT_NEAR(ewBefore - 3.14, G.totalEdgeWeight(), epsilon);
    } else {
        ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    ASSERT_EQ(1u, G.numberOfEdges());
    ASSERT_TRUE(G.hasEdge(0, 0));
    ASSERT_FALSE(G.hasEdge(0, 1));
    ASSERT_FALSE(G.hasEdge(2, 1));

    // test remove self-loop
    G.addEdge(2, 1);

    ewBefore = G.totalEdgeWeight();
    G.removeEdge(0, 0);
    if (G.isWeighted()) {
        ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    } else {
        ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    ASSERT_EQ(1u, G.numberOfEdges());
    ASSERT_FALSE(G.hasEdge(0, 0));
    ASSERT_FALSE(G.hasEdge(0, 1));
    ASSERT_TRUE(G.hasEdge(2, 1));

    // test from removeselfloops adapted for removeEdge
    G = createGraph(2);

    ewBefore = G.totalEdgeWeight();

    G.addEdge(0, 1);
    G.addEdge(0, 0, 3.14);
    G.addEdge(1, 1);

    if (G.isWeighted()) {
        EXPECT_NEAR(ewBefore + 3.14 + 2 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    } else {
        EXPECT_NEAR(ewBefore + 3 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    EXPECT_EQ(3u, G.numberOfEdges());
    EXPECT_TRUE(G.hasEdge(0, 0));
    EXPECT_TRUE(G.hasEdge(0, 1));
    EXPECT_TRUE(G.hasEdge(1, 1));
    EXPECT_EQ(G.numberOfSelfLoops(), 2u);

    // remove self-loops
    ewBefore = G.totalEdgeWeight();

    G.removeEdge(0, 0);
    G.removeEdge(1, 1);

    if (G.isWeighted()) {
        EXPECT_NEAR(ewBefore - defaultEdgeWeight - 3.14, G.totalEdgeWeight(), epsilon);
    } else {
        EXPECT_NEAR(ewBefore - 2 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon)
            << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
    }

    EXPECT_EQ(1u, G.numberOfEdges());
    EXPECT_FALSE(G.hasEdge(0, 0));
    EXPECT_FALSE(G.hasEdge(1, 1));
    EXPECT_TRUE(G.hasEdge(0, 1));
    EXPECT_EQ(0u, G.numberOfSelfLoops())
        << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
}

TEST_P(DHBGraphGTest, testRemoveAllEdges) {
    constexpr node graph_size = 100;
    DHBGraph G = createGraph(graph_size);
    for (size_t i = 0; i < graph_size; ++i) {
        for (size_t j = i + 1; j < graph_size; ++j) {
            G.addEdge(i, j);
        }
    }
    G.removeAllEdges();
    ASSERT_EQ(0u, G.numberOfEdges());
}

TEST_P(DHBGraphGTest, testRemoveAdjacentEdges_outEdge) {
    DHBGraph G = this->Ghouse;
    auto condition = [&](node) { return true; };
    auto [removedEdges, removedSelfLoops] = G.removeAdjacentEdges(2, condition, false);
    ASSERT_FALSE(G.hasEdge(2, 1));
    ASSERT_FALSE(G.hasEdge(2, 4));
    ASSERT_EQ(2, removedEdges);
    ASSERT_EQ(0, removedSelfLoops);
}

TEST_P(DHBGraphGTest, testRemoveAdjacentEdges_inEdge) {
    DHBGraph G = this->Ghouse;
    auto condition = [&](node) { return true; };

    auto [removedEdges, removedSelfLoops] = G.removeAdjacentEdges(2, condition, true);
    ASSERT_FALSE(G.hasEdge(0, 2));
    ASSERT_FALSE(G.hasEdge(3, 2));
    ASSERT_EQ(2, removedEdges);
    ASSERT_EQ(0, removedSelfLoops);
}

TEST_P(DHBGraphGTest, testRemoveAdjacentEdges_selfloop_outEdge) {
    DHBGraph G = this->Ghouse;
    G.addEdge(0, 0);
    ASSERT_EQ(1, G.numberOfSelfLoops());

    auto condition = [&](node) { return true; };
    auto [removedEdges, removedSelfLoops] = G.removeAdjacentEdges(0, condition, false);
    ASSERT_FALSE(G.hasEdge(0, 2));
    ASSERT_FALSE(G.hasEdge(0, 0));
    // Now, self loop is not considered in counting removed edges
    ASSERT_EQ(1, removedEdges);
    ASSERT_EQ(1, removedSelfLoops);
    ASSERT_EQ(0, G.numberOfSelfLoops());
}

TEST_P(DHBGraphGTest, testRemoveSelfLoops) {
    Aux::setNumberOfThreads(1);
    DHBGraph G = this->Ghouse;
    G.addEdge(0, 0);
    G.addEdge(1, 1);
    G.addEdge(2, 2);
    ASSERT_EQ(3, G.numberOfSelfLoops());
    G.removeSelfLoops();
    ASSERT_EQ(0, G.numberOfSelfLoops());
    ASSERT_FALSE(G.hasEdge(0, 0));
    ASSERT_FALSE(G.hasEdge(1, 1));
    ASSERT_FALSE(G.hasEdge(2, 2));
}

TEST_P(DHBGraphGTest, testHasEdge) {
    auto containsEdge = [&](std::pair<node, node> e) {
        auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
        return it != this->houseEdgesOut.end();
    };

    for (node u = 0; u < this->Ghouse.upperNodeIdBound(); u++) {
        for (node v = 0; v < this->Ghouse.upperNodeIdBound(); v++) {
            auto edge = std::make_pair(u, v);
            auto edgeReverse = std::make_pair(v, u);
            bool hasEdge = containsEdge(edge);
            bool hasEdgeReverse = containsEdge(edgeReverse);
            if (this->Ghouse.isDirected()) {
                ASSERT_EQ(hasEdge, this->Ghouse.hasEdge(u, v));
            } else {
                ASSERT_EQ(hasEdge || hasEdgeReverse, this->Ghouse.hasEdge(u, v));
            }
        }
    }
}

TEST_P(DHBGraphGTest, testNumberOfNodes) {
    ASSERT_EQ(this->n_house, this->Ghouse.numberOfNodes());

    DHBGraph G1 = createGraph(0);
    ASSERT_EQ(0u, G1.numberOfNodes());
    G1.addNode();
    ASSERT_EQ(1u, G1.numberOfNodes());
    G1.addNode();
    ASSERT_EQ(2u, G1.numberOfNodes());
}

TEST_P(DHBGraphGTest, testNumberOfEdges) {
    ASSERT_EQ(this->m_house, this->Ghouse.numberOfEdges());

    DHBGraph G1 = createGraph(5);
    ASSERT_EQ(0u, G1.numberOfEdges());
    G1.addEdge(0, 1);
    ASSERT_EQ(1u, G1.numberOfEdges());
    G1.addEdge(1, 2);
    ASSERT_EQ(2u, G1.numberOfEdges());
    G1.removeEdge(0, 1);
    ASSERT_EQ(1u, G1.numberOfEdges());
    G1.removeEdge(1, 2);
    ASSERT_EQ(0u, G1.numberOfEdges());
}

TEST_P(DHBGraphGTest, testNumberOfSelfLoops) {
    DHBGraph G = createGraph(3);
    G.addEdge(0, 1);
    ASSERT_EQ(0u, G.numberOfSelfLoops());
    G.addEdge(0, 0);
    ASSERT_EQ(1u, G.numberOfSelfLoops());
    G.addEdge(1, 1);
    G.addEdge(1, 2);
    ASSERT_EQ(2u, G.numberOfSelfLoops());
    G.removeEdge(0, 0);
    ASSERT_EQ(1u, G.numberOfSelfLoops());
}

TEST_P(DHBGraphGTest, DISABLED_testSelfLoopConversion) {
    Aux::Random::setSeed(1, false);
    count const runs = 100;
    const count n_max = 200;
    for (NetworKit::index i = 0; i < runs; i++) {
        bool directed = Aux::Random::probability() < 0.5;
        count n = Aux::Random::integer(n_max);
        DHBGraph G(n, false, directed);

        G.forNodes([&](node v) {
            double p = Aux::Random::probability();

            if (p < 0.1) { // new node
                n++;
                G.addNode();
            } else {                                     // new edge
                node u = Aux::Random::integer(v, n - 1); // self-loops possible
                G.addEdge(v, u);
            }
        });
        count measuredSelfLoops = countSelfLoopsManually(G);
        EXPECT_EQ(G.numberOfSelfLoops(), measuredSelfLoops);
        DHBGraph G_converted(G, false, !directed);
        EXPECT_EQ(G_converted.numberOfSelfLoops(), measuredSelfLoops);
    }
}

TEST_P(DHBGraphGTest, testUpperNodeIdBound) {
    ASSERT_EQ(5u, this->Ghouse.upperNodeIdBound());

    DHBGraph G1 = createGraph(0);
    ASSERT_EQ(0u, G1.upperNodeIdBound());
    G1.addNode();
    ASSERT_EQ(1u, G1.upperNodeIdBound());
    G1.addNode();
    ASSERT_EQ(2u, G1.upperNodeIdBound());
    G1.removeNode(1);
    ASSERT_EQ(2u, G1.upperNodeIdBound());
    G1.addNode();
    ASSERT_EQ(3u, G1.upperNodeIdBound());
}

/** EDGE ATTRIBUTES **/

TEST_P(DHBGraphGTest, testWeight) {
    this->Ghouse.forNodes([&](node u) {
        this->Ghouse.forNodes(
            [&](node v) { ASSERT_EQ(this->Ahouse[u][v], this->Ghouse.weight(u, v)); });
    });
}

TEST_P(DHBGraphGTest, testSetWeight) {
    DHBGraph G = createGraph(10);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    if (isWeighted()) {
        // edges should get weight defaultWeight on creation and setWeight should
        // overwrite this
        G.setWeight(1, 2, 2.718);
        EXPECT_EQ(defaultEdgeWeight, G.weight(0, 1));
        EXPECT_EQ(2.718, G.weight(1, 2));
        if (isDirected()) {
            EXPECT_EQ(nullWeight, G.weight(1, 0));
            EXPECT_EQ(nullWeight, G.weight(2, 1));
        } else {
            // undirected graph is symmetric
            EXPECT_EQ(defaultEdgeWeight, G.weight(1, 0));
            EXPECT_EQ(2.718, G.weight(2, 1));
        }

        // setting an edge weight should create the edge if it doesn't exists
        ASSERT_FALSE(G.hasEdge(5, 6));
        G.setWeight(5, 6, 56.0);
        ASSERT_EQ(56.0, G.weight(5, 6));
        ASSERT_EQ(isDirected() ? nullWeight : 56.0, G.weight(6, 5));
        ASSERT_TRUE(G.hasEdge(5, 6));

        // directed graphs are not symmetric, undirected are
        G.setWeight(2, 1, 5.243);
        if (isDirected()) {
            EXPECT_EQ(2.718, G.weight(1, 2));
            EXPECT_EQ(5.243, G.weight(2, 1));
        } else {
            EXPECT_EQ(5.243, G.weight(1, 2));
            EXPECT_EQ(5.243, G.weight(2, 1));
        }

        // self-loop
        G.addEdge(4, 4, 2.5);
        ASSERT_EQ(2.5, G.weight(4, 4));
        G.setWeight(4, 4, 3.14);
        ASSERT_EQ(3.14, G.weight(4, 4));
    }
}

TEST_P(DHBGraphGTest, testIncreaseWeight) {
    DHBGraph G = createGraph(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(3, 4, 3.14);

    if (G.isWeighted()) {
        G.increaseWeight(1, 2, 0.5);
        G.increaseWeight(3, 4, -0.5);

        ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
        ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(1, 2));
        ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));

        G.increaseWeight(2, 3, 0.3);
        ASSERT_EQ(0.3, G.weight(2, 3));
        ASSERT_TRUE(G.hasEdge(2, 3));

        if (G.isDirected()) {
            // reverse edges do net exist => weight should be nullWeight
            ASSERT_EQ(nullWeight, G.weight(1, 0));
            ASSERT_EQ(nullWeight, G.weight(2, 1));
            ASSERT_EQ(nullWeight, G.weight(4, 3));
        } else {
            ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
            ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(2, 1));
            ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));
        }
    } else {
        EXPECT_ANY_THROW(G.increaseWeight(1, 2, 0.3));
    }
}

/** SUMS **/

TEST_P(DHBGraphGTest, testTotalEdgeWeight) {
    DHBGraph G1 = createGraph(5);
    DHBGraph G2 = createGraph(5);
    G2.addEdge(0, 1, 3.14);

    if (this->Ghouse.isWeighted()) {
        ASSERT_EQ(0.0, G1.totalEdgeWeight());
        ASSERT_EQ(3.14, G2.totalEdgeWeight());
        ASSERT_EQ(36.0, this->Ghouse.totalEdgeWeight());
    } else {
        ASSERT_EQ(0 * defaultEdgeWeight, G1.totalEdgeWeight());
        ASSERT_EQ(1 * defaultEdgeWeight, G2.totalEdgeWeight());
        ASSERT_EQ(8 * defaultEdgeWeight, this->Ghouse.totalEdgeWeight());
    }
}

/** Collections **/

TEST_P(DHBGraphGTest, testNodeIterator) {
    auto testForward = [](const DHBGraph &G) {
        auto preIter = G.nodeRange().begin();
        auto postIter = G.nodeRange().begin();

        G.forNodes([&](const node u) {
            ASSERT_EQ(*preIter, u);
            ASSERT_EQ(*postIter, u);
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(preIter, G.nodeRange().end());
        ASSERT_EQ(postIter, G.nodeRange().end());

        DHBGraph G1(G);

        for (const auto u : DHBGraph::NodeRange(G)) {
            ASSERT_TRUE(G1.hasNode(u));
            G1.removeNode(u);
        }

        G1.forNodes([&](node u) { ASSERT_EQ(G1.degree(u), 0); });
    };

    auto testBackward = [](const DHBGraph &G) {
        const std::vector<node> nodes(DHBGraph::NodeRange(G).begin(), DHBGraph::NodeRange(G).end());
        std::vector<node> v;
        G.forNodes([&](node u) { v.push_back(u); });

        ASSERT_EQ(std::unordered_set<node>(nodes.begin(), nodes.end()).size(), nodes.size());
        ASSERT_EQ(nodes.size(), G.numberOfNodes());

        auto preIter = G.nodeRange().begin();
        auto postIter = G.nodeRange().begin();
        for (count i = 0; i < G.numberOfNodes(); ++i) {
            ++preIter;
            postIter++;
        }

        ASSERT_EQ(preIter, G.nodeRange().end());
        ASSERT_EQ(postIter, G.nodeRange().end());
        auto vecIter = nodes.rbegin();
        while (vecIter != nodes.rend()) {
            ASSERT_EQ(*vecIter, *(--preIter));
            if (postIter != G.nodeRange().end()) {
                ASSERT_NE(*vecIter, *(postIter--));
            } else {
                postIter--;
            }
            ASSERT_EQ(*vecIter, *postIter);
            ++vecIter;
        }

        ASSERT_EQ(preIter, G.nodeRange().begin());
        ASSERT_EQ(postIter, G.nodeRange().begin());
    };

    testForward(Ghouse);
    testBackward(Ghouse);
}

TEST_P(DHBGraphGTest, testEdgeIterator) {
    auto testForward = [&](const DHBGraph &G) {
        DHBGraph G1(G);
        auto preIter = G.edgeRange().begin();
        auto postIter = G.edgeRange().begin();

        G.forEdges([&](node, node) {
            ASSERT_EQ(preIter, postIter);
            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            G1.removeEdge(edge.u, edge.v);
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
        ASSERT_EQ(preIter, G.edgeRange().end());
        ASSERT_EQ(postIter, G.edgeRange().end());

        G1 = G;
        for (const auto edge : DHBGraph::EdgeRange(G)) {
            ASSERT_TRUE(G1.hasEdge(edge.u, edge.v));
            G1.removeEdge(edge.u, edge.v);
        }

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto testForwardWeighted = [&](const DHBGraph &G) {
        DHBGraph G1(G);
        auto preIter = G.edgeWeightRange().begin();
        auto postIter = preIter;

        G.forEdges([&](node, node) {
            ASSERT_EQ(preIter, postIter);

            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            ASSERT_DOUBLE_EQ(G.weight(edge.u, edge.v), edge.weight);
            G1.removeEdge(edge.u, edge.v);
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
        ASSERT_EQ(preIter, G.edgeWeightRange().end());
        ASSERT_EQ(postIter, G.edgeWeightRange().end());

        G1 = G;
        for (const auto edge : DHBGraph::EdgeWeightRange(G)) {
            ASSERT_TRUE(G1.hasEdge(edge.u, edge.v));
            ASSERT_DOUBLE_EQ(G1.weight(edge.u, edge.v), edge.weight);
            G1.removeEdge(edge.u, edge.v);
        }

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto testBackward = [&](const DHBGraph &G) {
        DHBGraph G1(G);
        auto preIter = G.edgeRange().begin();
        auto postIter = preIter;
        G.forEdges([&](node, node) {
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(preIter, G.edgeRange().end());
        ASSERT_EQ(postIter, G.edgeRange().end());

        G.forEdges([&](node, node) {
            --preIter;
            postIter--;
            ASSERT_EQ(preIter, postIter);
            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            G1.removeEdge(edge.u, edge.v);
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto testBackwardWeighted = [&](const DHBGraph &G) {
        DHBGraph G1(G);
        auto preIter = G.edgeWeightRange().begin();
        auto postIter = preIter;
        G.forEdges([&](node, node) {
            ++preIter;
            postIter++;
        });

        G.forEdges([&](node, node) {
            --preIter;
            postIter--;
            ASSERT_EQ(preIter, postIter);

            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            ASSERT_DOUBLE_EQ(G.weight(edge.u, edge.v), edge.weight);
            G1.removeEdge(edge.u, edge.v);
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto doTests = [&](const DHBGraph &G) {
        testForward(G);
        testForwardWeighted(G);
        testBackward(G);
        testBackwardWeighted(G);
    };

    doTests(Ghouse);
}

TEST_P(DHBGraphGTest, testNeighborsIterators) {
    auto iter = this->Ghouse.neighborRange(1).begin();
    this->Ghouse.forNeighborsOf(1, [&](node v) {
        ASSERT_TRUE(*iter == v);
        ++iter;
    });
    ASSERT_TRUE(iter == this->Ghouse.neighborRange(1).end());

    if (this->Ghouse.isWeighted()) {
        auto iterW = this->Ghouse.weightNeighborRange(1).begin();
        this->Ghouse.forNeighborsOf(1, [&](node v, edgeweight w) {
            ASSERT_TRUE((*iterW).first == v);
            ASSERT_TRUE((*iterW).second == w);
            ++iterW;
        });
        ASSERT_TRUE(iterW == this->Ghouse.weightNeighborRange(1).end());
    }
}

TEST_P(DHBGraphGTest, testNeighborsIterator_pre_increment) {
    auto iter = this->Ghouse.neighborRange(1).begin();
    std::vector<node> nodes;
    this->Ghouse.forNeighborsOf(1, [&](node v) { nodes.push_back(v); });
    auto orig_iter = iter;
    ++iter;
    ASSERT_NE(iter, orig_iter);
    ASSERT_EQ(*iter, nodes[1]);
}

TEST_P(DHBGraphGTest, testNeighborsIterator_post_increment) {
    auto iter = this->Ghouse.neighborRange(1).begin();
    std::vector<node> nodes;
    this->Ghouse.forNeighborsOf(1, [&](node v) { nodes.push_back(v); });
    auto orig_iter = iter;
    auto temp = orig_iter++;
    ASSERT_EQ(*temp, nodes[0]);
    ASSERT_EQ(orig_iter, ++iter);
}

TEST_P(DHBGraphGTest, testNeighborsIterator_pre_decrement) {
    std::vector<node> nodes;
    this->Ghouse.forNeighborsOf(1, [&](node v) { nodes.push_back(v); });
    auto iter = this->Ghouse.neighborRange(1).begin();
    auto orig_iter = iter;
    ++iter;
    --iter;
    ASSERT_EQ(iter, orig_iter);
    ASSERT_EQ(*iter, nodes[0]);
}

TEST_P(DHBGraphGTest, testNeighborsIterator_post_decrement) {
    std::vector<node> nodes;
    this->Ghouse.forNeighborsOf(1, [&](node v) { nodes.push_back(v); });
    auto iter_post = this->Ghouse.neighborRange(1).begin();
    auto postfix_orig_iter = this->Ghouse.neighborRange(1).begin();
    ++iter_post;
    auto temp = iter_post--;
    ASSERT_EQ(*temp, nodes[1]);
    ASSERT_EQ(postfix_orig_iter, iter_post);
    ASSERT_EQ(*postfix_orig_iter, nodes[0]);
}

TEST_P(DHBGraphGTest, testGetIthNeighbor) {
    DHBGraph G = createGraph(4);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(0, 3);

    node neighbor_0 = G.getIthNeighbor(0, 0);
    node neighbor_1 = G.getIthNeighbor(0, 1);
    node neighbor_2 = G.getIthNeighbor(0, 2);
    node neighbor_3 = G.getIthNeighbor(0, 3);

    ASSERT_EQ(1, neighbor_0);
    ASSERT_EQ(2, neighbor_1);
    ASSERT_EQ(3, neighbor_2);
    ASSERT_EQ(none, neighbor_3);
}

TEST_P(DHBGraphGTest, testGetIthNeighborWeight) {
    DHBGraph G = createGraph(4);
    G.addEdge(0, 1, 2.f);
    G.addEdge(0, 2, 2.f);
    G.addEdge(0, 3, 2.f);

    node neighbor_0_weight = G.getIthNeighborWeight(0, 0);
    node neighbor_1_weight = G.getIthNeighborWeight(0, 1);
    node neighbor_2_weight = G.getIthNeighborWeight(0, 2);
    node neighbor_3_weight = G.getIthNeighborWeight(0, 3);

    ASSERT_EQ(G.isWeighted() ? 2.f : defaultEdgeWeight, neighbor_0_weight);
    ASSERT_EQ(G.isWeighted() ? 2.f : defaultEdgeWeight, neighbor_1_weight);
    ASSERT_EQ(G.isWeighted() ? 2.f : defaultEdgeWeight, neighbor_2_weight);
    ASSERT_EQ(defaultEdgeWeight, neighbor_3_weight);
}

TEST_P(DHBGraphGTest, testGetIthNeighborWithWeight) {
    DHBGraph G = createGraph(4);
    G.addEdge(0, 1, 2.f);
    G.addEdge(0, 2, 2.f);
    G.addEdge(0, 3, 2.f);

    auto [neighbor_0, weight_0] = G.getIthNeighborWithWeight(0, 0);
    auto [neighbor_1, weight_1] = G.getIthNeighborWithWeight(0, 1);
    auto [neighbor_2, weight_2] = G.getIthNeighborWithWeight(0, 2);
    auto [neighbor_3, weight_3] = G.getIthNeighborWithWeight(0, 3);

    ASSERT_EQ(1, neighbor_0);
    ASSERT_EQ(2, neighbor_1);
    ASSERT_EQ(3, neighbor_2);
    ASSERT_EQ(none, neighbor_3);

    ASSERT_EQ(G.isWeighted() ? 2.f : defaultEdgeWeight, weight_0);
    ASSERT_EQ(G.isWeighted() ? 2.f : defaultEdgeWeight, weight_1);
    ASSERT_EQ(G.isWeighted() ? 2.f : defaultEdgeWeight, weight_2);
    ASSERT_EQ(defaultEdgeWeight, weight_3);
}

TEST_P(DHBGraphGTest, testGetIthNeighborWithId) {
    DHBGraph G = createGraph(4);
    G.addEdge(0, 1, 2.f);
    G.addEdge(0, 2, 2.f);
    G.addEdge(0, 3, 2.f);
    G.indexEdges();

    auto [neighbor_0, id_0] = G.getIthNeighborWithId(0, 0);
    auto [neighbor_1, id_1] = G.getIthNeighborWithId(0, 1);
    auto [neighbor_2, id_2] = G.getIthNeighborWithId(0, 2);
    auto [neighbor_3, id_3] = G.getIthNeighborWithId(0, 3);

    ASSERT_EQ(1, neighbor_0);
    ASSERT_EQ(2, neighbor_1);
    ASSERT_EQ(3, neighbor_2);

    ASSERT_EQ(0, id_0);
    ASSERT_EQ(1, id_1);
    ASSERT_EQ(2, id_2);
    ASSERT_EQ(none, id_3);
}

TEST_P(DHBGraphGTest, testSetWeightAtIthNeighbor) {
    DHBGraph G = createGraph(4);
    G.addEdge(0, 1, 2.f);
    G.addEdge(0, 2, 2.f);
    G.addEdge(0, 3, 2.f);
    G.setWeightAtIthNeighbor(unsafe, 0, 0, 5.f);
    double const result = G.isWeighted() ? 5.f : defaultEdgeWeight;
    ASSERT_EQ(result, G.getIthNeighborWeight(0, 0));
}

TEST_P(DHBGraphGTest, testIndexOfNeighbor) {
    DHBGraph G = createGraph(5);
    G.addEdge(0, 1, 2.f);
    G.addEdge(0, 2, 2.f);
    G.addEdge(0, 3, 2.f);

    NetworKit::index i_1 = G.indexOfNeighbor(0, 2);
    NetworKit::index i_2 = G.indexOfNeighbor(0, 3);
    NetworKit::index i_none = G.indexOfNeighbor(0, 4);

    ASSERT_EQ(1, i_1);
    ASSERT_EQ(2, i_2);
    ASSERT_EQ(none, i_none);
}

TEST_P(DHBGraphGTest, testIndexOfNeighbor_2) {
    // Create a larger graph with more nodes
    DHBGraph G = createGraph(10); // Adjust the number of nodes as needed

    // Add multiple edges
    G.addEdge(0, 1, 2.f);
    G.addEdge(0, 2, 2.f);
    G.addEdge(0, 3, 2.f);
    G.addEdge(0, 4, 2.f);
    G.addEdge(0, 5, 2.f);
    G.addEdge(0, 6, 2.f);
    G.addEdge(0, 7, 2.f);
    G.addEdge(0, 8, 2.f);
    G.addEdge(0, 9, 2.f);

    // Test index retrieval for existing neighbors
    NetworKit::index i_1 = G.indexOfNeighbor(0, 2);
    NetworKit::index i_2 = G.indexOfNeighbor(0, 4);
    NetworKit::index i_3 = G.indexOfNeighbor(0, 6);
    NetworKit::index i_4 = G.indexOfNeighbor(0, 9);

    // Test index retrieval for non-existing neighbor
    NetworKit::index i_none = G.indexOfNeighbor(0, 10); // Node 10 does not exist

    // Check expected indices
    ASSERT_EQ(1, i_1);
    ASSERT_EQ(3, i_2);
    ASSERT_EQ(5, i_3);
    ASSERT_EQ(8, i_4);
    ASSERT_EQ(none, i_none);
}

/** NODE ITERATORS **/

TEST_P(DHBGraphGTest, testForNodes) {
    DHBGraph G = createGraph(3);
    std::vector<bool> visited(4, false);
    G.forNodes([&](node v) {
        ASSERT_FALSE(visited[v]);
        if (v == 2) {
            G.addNode();
        }
        visited[v] = true;
    });
    for (bool b : visited) {
        ASSERT_TRUE(b);
    }
}

TEST_P(DHBGraphGTest, testParallelForNodes) {
    std::vector<node> visited(Ghouse.upperNodeIdBound());
    this->Ghouse.parallelForNodes([&](node u) { visited[u] = u; });

    Aux::Parallel::sort(visited.begin(), visited.end());

    ASSERT_EQ(5u, visited.size());
    for (NetworKit::index i = 0; i < this->Ghouse.upperNodeIdBound(); i++) {
        ASSERT_EQ(i, visited[i]);
    }
}

TEST_P(DHBGraphGTest, forNodesWhile) {
    count n = 100;
    DHBGraph G = createGraph(n);
    count stopAfter = 10;
    count nodesSeen = 0;

    G.forNodesWhile([&]() { return nodesSeen < stopAfter; }, [&](node) { nodesSeen++; });

    ASSERT_EQ(stopAfter, nodesSeen);
}

TEST_P(DHBGraphGTest, testForNodesInRandomOrder) {
    count n = 1000;
    count samples = 100;
    double maxAbsoluteError = 0.005;
    DHBGraph G = createGraph(n);

    node lastNode = n / 2;
    count greaterLastNode = 0;
    count smallerLastNode = 0;
    std::vector<count> visitCount(n, 0);

    for (count i = 0; i < samples; i++) {
        G.forNodesInRandomOrder([&](node v) {
            if (v > lastNode) {
                greaterLastNode++;
            } else {
                smallerLastNode++;
            }
            visitCount[v]++;
            lastNode = v;
        });
    }

    for (node v = 0; v < n; v++) {
        ASSERT_EQ(samples, visitCount[v]);
    }

    ASSERT_NEAR(0.5, (double)greaterLastNode / n / samples, maxAbsoluteError);
    ASSERT_NEAR(0.5, (double)smallerLastNode / n / samples, maxAbsoluteError);
}

TEST_P(DHBGraphGTest, testForNodePairs) {
    count n = 10;
    count m = n * (n - 1) / 2;
    DHBGraph G = createGraph(n);

    // add all edges
    G.forNodePairs([&](node u, node v) {
        ASSERT_FALSE(G.hasEdge(u, v));
        G.addEdge(u, v);
        ASSERT_TRUE(G.hasEdge(u, v));
    });

    EXPECT_EQ(m, G.numberOfEdges());

    // remove all edges
    G.forNodePairs([&](node u, node v) {
        ASSERT_TRUE(G.hasEdge(u, v));
        G.removeEdge(u, v);
        ASSERT_FALSE(G.hasEdge(u, v));
    });

    EXPECT_EQ(0u, G.numberOfEdges());
}

TEST_P(DHBGraphGTest, testForNodePairsParallel) {
    count n = 10;
    DHBGraph G = createGraph(n);

    G.forNodePairsParallel([&](node u, node v) { ASSERT_FALSE(G.hasEdge(u, v)); });

    for (node u = 0; u < n; ++u) {
        for (node v = u + 1; v < n; ++v) {
            G.addEdge(u, v);
        }
    }

    G.forNodePairsParallel([&](node u, node v) { ASSERT_TRUE(G.hasEdge(u, v)); });

    count m = n * (n - 1) / 2;
    EXPECT_EQ(m, G.numberOfEdges());
}

/** EDGE ITERATORS **/

TEST_P(DHBGraphGTest, testForEdges) {
    DHBGraph G = createGraph(4);
    G.addEdge(0, 1); // 0 * 1 = 0
    G.addEdge(1, 2); // 1 * 2 = 2
    G.addEdge(3, 2); // 3 * 2 = 1 (mod 5)
    G.addEdge(2, 2); // 2 * 2 = 4
    G.addEdge(3, 1); // 3 * 1 = 3

    std::vector<bool> edgesSeen(5, false);

    G.forEdges([&](node u, node v) {
        ASSERT_TRUE(G.hasEdge(u, v));
        NetworKit::index id = (u * v) % 5;
        edgesSeen[id] = true;
    });

    for (auto b : edgesSeen) {
        ASSERT_TRUE(b);
    }
}

TEST_P(DHBGraphGTest, testForWeightedEdges) {
    double epsilon = 1e-6;

    DHBGraph G = createGraph(4);
    G.addEdge(0, 1, 0.1); // 0 * 1 = 0
    G.addEdge(3, 2, 0.2); // 3 * 2 = 1 (mod 5)
    G.addEdge(1, 2, 0.3); // 1 * 2 = 2
    G.addEdge(3, 1, 0.4); // 3 * 1 = 3
    G.addEdge(2, 2, 0.5); // 2 * 2 = 4

    std::vector<bool> edgesSeen(5, false);

    edgeweight weightSum = 0;
    G.forEdges([&](node u, node v, edgeweight ew) {
        ASSERT_TRUE(G.hasEdge(u, v));
        ASSERT_EQ(G.weight(u, v), ew);

        NetworKit::index id = (u * v) % 5;
        edgesSeen[id] = true;
        if (G.isWeighted()) {
            ASSERT_NEAR((id + 1) * 0.1, ew, epsilon);
        } else {
            ASSERT_EQ(defaultEdgeWeight, ew);
        }
        weightSum += ew;
    });

    for (auto b : edgesSeen) {
        ASSERT_TRUE(b);
    }
    if (G.isWeighted()) {
        ASSERT_NEAR(1.5, weightSum, epsilon);
    } else {
        ASSERT_NEAR(5 * defaultEdgeWeight, weightSum, epsilon);
    }
}

TEST_P(DHBGraphGTest, testParallelForWeightedEdges) {
    count n = 4;
    DHBGraph G = createGraph(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v, 1.0); });

    edgeweight weightSum = 0.0;
    G.parallelForEdges([&](node, node, edgeweight ew) {
#pragma omp atomic
        weightSum += ew;
    });

    ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
}

TEST_P(DHBGraphGTest, testParallelForEdges) {
    count n = 4;
    DHBGraph G = createGraph(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v); });

    edgeweight weightSum = 0.0;
    G.parallelForEdges([&](node, node) {
#pragma omp atomic
        weightSum += 1;
    });

    ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
}

/** NEIGHBORHOOD ITERATORS **/

TEST_P(DHBGraphGTest, testForNeighborsOf) {
    std::vector<std::pair<node, edgeweight>> visited;
    this->Ghouse.forNeighborsOf(
        3, [&](node u, edgeweight ew) { visited.push_back(std::make_pair(u, ew)); });

    // should sort after the first element
    Aux::Parallel::sort(visited.begin(), visited.end());

    if (isGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        ASSERT_EQ(defaultEdgeWeight, visited[1].second);
        ASSERT_EQ(defaultEdgeWeight, visited[2].second);
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(1.0, visited[0].second);
        ASSERT_EQ(7.0, visited[1].second);
        ASSERT_EQ(6.0, visited[2].second);
    }

    if (isDirectedGraph()) {
        ASSERT_EQ(2u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        ASSERT_EQ(defaultEdgeWeight, visited[1].second);
    }

    if (isWeightedDirectedGraph()) {
        ASSERT_EQ(2u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(1.0, visited[0].second);
        ASSERT_EQ(7.0, visited[1].second);
    }
}

TEST_P(DHBGraphGTest, testForNeighborsOfParallel) {
    std::vector<std::pair<node, edgeweight>> visited;
    std::mutex visitedMutex;

    this->Ghouse.forNeighborsOfParallel(3, [&](node u, edgeweight ew) {
        std::lock_guard<std::mutex> lock(visitedMutex);
        visited.push_back(std::make_pair(u, ew));
    });

    // should sort after the first element
    std::sort(visited.begin(), visited.end());

    if (isGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        ASSERT_EQ(defaultEdgeWeight, visited[1].second);
        ASSERT_EQ(defaultEdgeWeight, visited[2].second);
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(1.0, visited[0].second);
        ASSERT_EQ(7.0, visited[1].second);
        ASSERT_EQ(6.0, visited[2].second);
    }

    if (isDirectedGraph()) {
        ASSERT_EQ(2u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        ASSERT_EQ(defaultEdgeWeight, visited[1].second);
    }

    if (isWeightedDirectedGraph()) {
        ASSERT_EQ(2u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(1.0, visited[0].second);
        ASSERT_EQ(7.0, visited[1].second);
    }
}

TEST_P(DHBGraphGTest, testForEdgesOf) {
    count m = 0;
    std::vector<int> visited(this->m_house, 0);

    this->Ghouse.forNodes([&](node u) {
        this->Ghouse.forEdgesOf(u, [&](node v, node w) {
            // edges should be v to w, so if we iterate over edges from u, u should be
            // equal v
            EXPECT_EQ(u, v);

            auto e = std::make_pair(v, w);
            auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
            if (!isDirected() && it == this->houseEdgesOut.end()) {
                auto e2 = std::make_pair(w, v);
                it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e2);
            }

            EXPECT_TRUE(it != this->houseEdgesOut.end());

            // find index in edge array
            int i = std::distance(this->houseEdgesOut.begin(), it);
            if (isDirected()) {
                // make sure edge was not visited before (would be visited twice)
                EXPECT_EQ(0, visited[i]);
            }

            // mark edge as visited
            visited[i]++;
            m++;
        });
    });

    if (isDirected()) {
        // we iterated over all outgoing edges once
        EXPECT_EQ(this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(1, c);
        }
    } else {
        // we iterated over all edges in both directions
        EXPECT_EQ(2 * this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(2, c);
        }
    }
}

TEST_P(DHBGraphGTest, testForWeightedEdgesOf) {
    count m = 0;
    std::vector<int> visited(this->m_house, 0);
    double sumOfWeights = 0;

    this->Ghouse.forNodes([&](node u) {
        this->Ghouse.forEdgesOf(u, [&](node v, node w, edgeweight ew) {
            // edges should be v to w, so if we iterate over edges from u, u should be
            // equal v
            EXPECT_EQ(u, v);
            sumOfWeights += ew;
            auto e = std::make_pair(v, w);
            auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
            if (!isDirected() && it == this->houseEdgesOut.end()) {
                auto e2 = std::make_pair(w, v);
                it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e2);
            }

            EXPECT_TRUE(it != this->houseEdgesOut.end());

            // find index in edge array
            int i = std::distance(this->houseEdgesOut.begin(), it);
            if (isDirected()) {
                // make sure edge was not visited before (would be visited twice)
                EXPECT_EQ(0, visited[i]);
            }

            // mark edge as visited
            visited[i]++;
            m++;
        });
    });

    if (isGraph()) {
        EXPECT_EQ(sumOfWeights, m);
        EXPECT_EQ(2 * this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(2, c);
        }
    }

    if (isWeightedGraph()) {
        // we iterated over all edges in both directions
        EXPECT_EQ(2 * this->m_house, m);
        EXPECT_EQ(sumOfWeights, 72);
        for (auto c : visited) {
            EXPECT_EQ(2, c);
        }
    }

    if (isDirectedGraph()) {
        // we iterated over all outgoing edges once
        EXPECT_EQ(this->m_house, m);
        EXPECT_EQ(sumOfWeights, m);
        for (auto c : visited) {
            EXPECT_EQ(1, c);
        }
    }

    if (isWeightedDirectedGraph()) {
        EXPECT_EQ(sumOfWeights, 36);
        EXPECT_EQ(this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(1, c);
        }
    }
}

TEST_P(DHBGraphGTest, testForInNeighborsOfParallel) {
    for (bool setIndexed : {false, true}) {
        if (setIndexed) {
            this->Ghouse.indexEdges();
        }
        std::vector<std::pair<node, edgeweight>> visited;
        std::mutex visitedMutex;

        this->Ghouse.forInNeighborsOfParallel(3, [&](node v, edgeweight ew) {
            std::lock_guard<std::mutex> lock(visitedMutex);
            visited.push_back({v, ew});
        });

        std::sort(visited.begin(), visited.end());

        if (isGraph()) {
            ASSERT_EQ(3u, visited.size());
            ASSERT_EQ(1u, visited[0].first);
            ASSERT_EQ(2u, visited[1].first);
            ASSERT_EQ(4u, visited[2].first);
            ASSERT_EQ(defaultEdgeWeight, visited[0].second);
            ASSERT_EQ(defaultEdgeWeight, visited[1].second);
            ASSERT_EQ(defaultEdgeWeight, visited[2].second);
        }

        if (isWeightedGraph()) {
            ASSERT_EQ(3u, visited.size());
            ASSERT_EQ(1u, visited[0].first);
            ASSERT_EQ(2u, visited[1].first);
            ASSERT_EQ(4u, visited[2].first);
            ASSERT_EQ(1.0, visited[0].second);
            ASSERT_EQ(7.0, visited[1].second);
            ASSERT_EQ(6.0, visited[2].second);
        }

        if (isDirectedGraph()) {
            ASSERT_EQ(1u, visited.size());
            ASSERT_EQ(4u, visited[0].first);
            ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        }

        if (isWeightedDirectedGraph()) {
            ASSERT_EQ(1u, visited.size());
            ASSERT_EQ(4u, visited[0].first);
            ASSERT_EQ(6.0, visited[0].second);
        }
    }
}

TEST_P(DHBGraphGTest, testForInNeighborsOf) {
    for (bool setIndexed : {false, true}) {
        if (setIndexed) {
            this->Ghouse.indexEdges();
        }

        std::vector<std::pair<node, edgeweight>> visited;
        this->Ghouse.forInNeighborsOf(3, [&](node v, edgeweight ew) {
            visited.push_back({v, ew});
        });
        std::sort(visited.begin(), visited.end());

        if (isGraph()) {
            ASSERT_EQ(3u, visited.size());
            ASSERT_EQ(1u, visited[0].first);
            ASSERT_EQ(2u, visited[1].first);
            ASSERT_EQ(4u, visited[2].first);
            ASSERT_EQ(defaultEdgeWeight, visited[0].second);
            ASSERT_EQ(defaultEdgeWeight, visited[1].second);
            ASSERT_EQ(defaultEdgeWeight, visited[2].second);
        }

        if (isWeightedGraph()) {
            ASSERT_EQ(3u, visited.size());
            ASSERT_EQ(1u, visited[0].first);
            ASSERT_EQ(2u, visited[1].first);
            ASSERT_EQ(4u, visited[2].first);
            ASSERT_EQ(1.0, visited[0].second);
            ASSERT_EQ(7.0, visited[1].second);
            ASSERT_EQ(6.0, visited[2].second);
        }

        if (isDirectedGraph()) {
            ASSERT_EQ(1u, visited.size());
            ASSERT_EQ(4u, visited[0].first);
            ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        }

        if (isWeightedDirectedGraph()) {
            ASSERT_EQ(1u, visited.size());
            ASSERT_EQ(4u, visited[0].first);
            ASSERT_EQ(6.0, visited[0].second);
        }
    }
}

TEST_P(DHBGraphGTest, testForInEdgesOf) {
    for (bool setIndexed : {false, true}) {
        if (setIndexed) {
            this->Ghouse.indexEdges();
        }

        // add self-loop
        this->Ghouse.addEdge(3, 3, 2.5);
        this->Ahouse[3][3] = 2.5;

        std::vector<edgeweight> visited(this->n_house, -1.0);
        this->Ghouse.forInEdgesOf(3, [&](node, node u, edgeweight ew) {
            ASSERT_EQ(-1.0, visited[u]);
            visited[u] = ew;
        });

        if (isGraph()) {
            ASSERT_EQ(-1.0, visited[0]);
            ASSERT_EQ(defaultEdgeWeight, visited[1]);
            ASSERT_EQ(defaultEdgeWeight, visited[2]);
            ASSERT_EQ(defaultEdgeWeight, visited[3]);
            ASSERT_EQ(defaultEdgeWeight, visited[4]);
        }

        if (isWeightedGraph()) {
            ASSERT_EQ(-1.0, visited[0]);
            ASSERT_EQ(this->Ahouse[3][1], visited[1]);
            ASSERT_EQ(this->Ahouse[3][2], visited[2]);
            ASSERT_EQ(this->Ahouse[3][3], visited[3]);
            ASSERT_EQ(this->Ahouse[3][4], visited[4]);
        }

        if (isDirectedGraph()) {
            ASSERT_EQ(-1.0, visited[0]);
            ASSERT_EQ(-1.0, visited[1]);
            ASSERT_EQ(-1.0, visited[2]);
            ASSERT_EQ(defaultEdgeWeight, visited[3]);
            ASSERT_EQ(defaultEdgeWeight, visited[4]);
        }

        if (isWeightedDirectedGraph()) {
            ASSERT_EQ(-1.0, visited[0]);
            ASSERT_EQ(-1.0, visited[1]);
            ASSERT_EQ(-1.0, visited[2]);
            ASSERT_EQ(this->Ahouse[3][3], visited[3]);
            ASSERT_EQ(this->Ahouse[4][3], visited[4]);
        }
    }
}

/** REDUCTION ITERATORS **/

TEST_P(DHBGraphGTest, testParallelSumForNodes) {
    count n = 10;
    DHBGraph G = createGraph(n);
    double sum = G.sumForNodesParallel([](node v) { return 2 * v + 0.5; });

    double expected_sum = n * (n - 1) + n * 0.5;
    ASSERT_EQ(expected_sum, sum);
}

TEST_P(DHBGraphGTest, testParallelSumForEdges) {
    for (bool setIndexed : {false, true}) {
        if (setIndexed) {
            this->Ghouse.indexEdges();
        }

        double sum =
            this->Ghouse.sumForEdgesParallel([](node, node, edgeweight ew) { return 1.5 * ew; });

        double expected_sum = 1.5 * this->Ghouse.totalEdgeWeight();
        ASSERT_EQ(expected_sum, sum);
    }
}

/** GRAPH SEARCHES **/

TEST_P(DHBGraphGTest, testEdgeIndexGenerationDirected) {
    DHBGraph G = DHBGraph(10, false, true);
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(2, 2);
    G.addEdge(5, 6);
    G.addEdge(6, 5);
    G.addEdge(1, 2);

    G.indexEdges();

    // Check consecutiveness of edgeids according to edge iterators
    edgeid expectedId = 0;
    G.forEdges([&](node u, node v) { EXPECT_EQ(expectedId++, G.edgeId(u, v)); });

    // Add some more edges
    EXPECT_EQ(6, G.upperEdgeIdBound());
    G.addEdge(8, 9);
    EXPECT_EQ(7, G.upperEdgeIdBound());
    G.addEdge(9, 8);

    // Check that asymmetric edges do not have the same id
    EXPECT_NE(G.edgeId(6, 5), G.edgeId(5, 6));
    EXPECT_NE(G.edgeId(2, 1), G.edgeId(1, 2));
    EXPECT_NE(G.edgeId(9, 8), G.edgeId(8, 9));
    EXPECT_EQ(7, G.edgeId(9, 8));
    EXPECT_EQ(8, G.upperEdgeIdBound());
}

TEST_P(DHBGraphGTest, testEdgeIndexGenerationUndirected) {
    DHBGraph G = DHBGraph(10, false, false);

    G.addEdge(0, 0);
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(2, 2);
    G.addEdge(5, 6);

    G.indexEdges();

    // Check consecutiveness of edgeids according to edge iterators
    edgeid expectedId = 0;
    G.forEdges([&](node u, node v) { EXPECT_EQ(expectedId++, G.edgeId(u, v)); });

    // Add some more edges. This will likely destroy consecutiveness...
    G.addEdge(3, 4);
    G.addEdge(7, 8);
    EXPECT_EQ(6, G.edgeId(7, 8));
    EXPECT_EQ(7, G.upperEdgeIdBound());

    // Anyway, check uniqueness and validity of the edgeids
    std::set<edgeid> ids;
    edgeid upperEdgeIdBound = G.upperEdgeIdBound();

    G.forEdges([&](node u, node v) {
        edgeid id = G.edgeId(u, v);
        EXPECT_EQ(id, G.edgeId(v, u));
        EXPECT_LT(id, upperEdgeIdBound);

        EXPECT_NE(none, id);
        EXPECT_FALSE(ids.erase(id));
        ids.insert(id);
    });
}

TEST_P(DHBGraphGTest, testEdgeIndexResolver) {
    DHBGraph G = createGraph(10);
    G.indexEdges();

    G.addEdge(0, 0);
    G.addEdge(5, 6);
    G.addEdge(2, 2);
    G.addEdge(4, 2);

    if (G.isDirected())
        G.addEdge(3, 2);

    std::map<std::pair<node, node>, edgeid> expectedEdges;
    expectedEdges[std::make_pair(0, 0)] = 0;
    expectedEdges[std::make_pair(5, 6)] = 1;
    expectedEdges[std::make_pair(2, 2)] = 2;
    expectedEdges[std::make_pair(4, 2)] = 3;
    expectedEdges[std::make_pair(2, 4)] = 3;
    expectedEdges[std::make_pair(3, 2)] = 4;

    G.forEdges([&](node, node, edgeid eid) {
        auto edge = G.edgeById(eid);
        EXPECT_EQ(expectedEdges[edge], eid);
    });
}

TEST_P(DHBGraphGTest, testForEdgesWithIds) {
    std::vector<DHBGraph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0);
        graph->addEdge(1, 2);
        graph->addEdge(4, 5);

        // No edge indices

        count m = 0;
        graph->forEdges([&](node, node, edgeid eid) {
            EXPECT_EQ(none, eid);
            m++;
        });
        ASSERT_EQ(3u, m);

        // With edge indices
        graph->indexEdges();

        edgeid expectedId = 0;
        m = 0;
        graph->forEdges([&](node, node, edgeid eid) {
            EXPECT_EQ(expectedId++, eid);
            EXPECT_LT(eid, graph->upperEdgeIdBound());
            m++;
        });
        ASSERT_EQ(3u, m);
    }
}

TEST_P(DHBGraphGTest, testForWeightedEdgesWithIds) {
    std::vector<DHBGraph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0, 2);
        graph->addEdge(1, 2, 2);
        graph->addEdge(4, 5, 2);

        // No edge indices

        count m = 0;
        edgeweight sum = 0;
        graph->forEdges([&](node, node, edgeweight ew, edgeid eid) {
            EXPECT_EQ(none, eid);
            m++;
            sum += ew;
        });
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }

        // With edge indices
        graph->indexEdges();

        edgeid expectedId = 0;
        m = 0;
        sum = .0;
        graph->forEdges([&](node, node, edgeweight ew, edgeid eid) {
            EXPECT_EQ(expectedId++, eid);
            EXPECT_LT(eid, graph->upperEdgeIdBound());
            m++;
            sum += ew;
        });
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }
    }
}

TEST_P(DHBGraphGTest, testParallelForEdgesWithIds) {
    std::vector<DHBGraph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0);
        graph->addEdge(1, 2);
        graph->addEdge(4, 5);

        // No edge indices
        count m = 0;
        edgeid sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeid eid) {
#pragma omp atomic
            m++;
            ASSERT_EQ(none, eid);
        });
        ASSERT_EQ(3u, m);

        // With edge indices
        graph->indexEdges();

        m = 0;
        edgeid expectedId = 0;
        sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeid eid) {
#pragma omp atomic
            expectedId++;
#pragma omp atomic
            sumedgeid += eid;
#pragma omp atomic
            m++;
        });
        ASSERT_EQ(expectedId, graph->upperEdgeIdBound());
        ASSERT_EQ(sumedgeid, ((graph->upperEdgeIdBound() - 1) * graph->upperEdgeIdBound()) / 2);
        ASSERT_EQ(3u, m);
    }
}

TEST_P(DHBGraphGTest, testParallelForWeightedEdgesWithIds) {
    std::vector<DHBGraph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0, 2);
        graph->addEdge(1, 2, 2);
        graph->addEdge(4, 5, 2);

        // No edge indices

        count m = 0;
        edgeweight sum = 0;
        edgeid sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeweight ew, edgeid eid) {
#pragma omp atomic
            m++;
#pragma omp atomic
            sum += ew;
            ASSERT_EQ(none, eid);
        });
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }

        // With edge indices
        graph->indexEdges();

        m = 0;
        sum = .0;
        edgeid expectedId = 0;
        sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeweight ew, edgeid eid) {
#pragma omp atomic
            expectedId++;
#pragma omp atomic
            sumedgeid += eid;
#pragma omp atomic
            m++;
#pragma omp atomic
            sum += ew;
        });
        ASSERT_EQ(expectedId, graph->upperEdgeIdBound());
        ASSERT_EQ(sumedgeid, ((graph->upperEdgeIdBound() - 1) * graph->upperEdgeIdBound()) / 2);
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }
    }
}

TEST_P(DHBGraphGTest, testSortEdges_default) {
    DHBGraph G = this->Ghouse;
    G.indexEdges();
    DHBGraph origG = G;
    G.sortEdges();
    std::vector<std::tuple<node, node, edgeweight, edgeid>> edges;
    std::vector<std::tuple<node, edgeweight, edgeid>> outEdges;
    edges.reserve(origG.numberOfEdges() * 4);
    origG.forNodes([&](node u) {
        origG.forEdgesOf(
            u, [&](node, node v, edgeweight w, edgeid eid) { outEdges.emplace_back(v, w, eid); });
        Aux::Parallel::sort(outEdges.begin(), outEdges.end());
        for (auto x : outEdges) {
            edges.emplace_back(u, std::get<0>(x), std::get<1>(x), std::get<2>(x));
        }
        outEdges.clear();
    });

    auto it = edges.begin();
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
            ASSERT_NE(it, edges.end());
            EXPECT_EQ(*it, std::make_tuple(u, v, w, eid));
            ++it;
        });
    });
}

TEST_P(DHBGraphGTest, testSortEdges_specified_lambda) {
    DHBGraph G = this->Ghouse;
    auto customSortLambda = [](node a_vertex, edgeweight, edgeid, node b_vertex, edgeweight,
                               edgeid) { return a_vertex < b_vertex; };

    G.indexEdges();
    DHBGraph origG = G;
    G.sortEdges(std::move(customSortLambda));
    std::vector<std::tuple<node, node, edgeweight, edgeid>> edges;
    std::vector<std::tuple<node, edgeweight, edgeid>> outEdges;
    edges.reserve(origG.numberOfEdges() * 4);
    origG.forNodes([&](node u) {
        origG.forEdgesOf(
            u, [&](node, node v, edgeweight w, edgeid eid) { outEdges.emplace_back(v, w, eid); });
        Aux::Parallel::sort(outEdges.begin(), outEdges.end());
        for (auto x : outEdges) {
            edges.emplace_back(u, std::get<0>(x), std::get<1>(x), std::get<2>(x));
        }
        outEdges.clear();
    });

    auto it = edges.begin();
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
            ASSERT_NE(it, edges.end());
            EXPECT_EQ(*it, std::make_tuple(u, v, w, eid));
            ++it;
        });
    });
}

TEST_P(DHBGraphGTest, testSwapEdges) {
    // For reference:
    // Edges of Ghouse (undirected, so we added the back edges in the second row):
    // {3, 1} {1, 3}
    // {1, 0} {0, 1}
    // {0, 2} {2, 0}
    // {2, 1} {1, 2}
    // {1, 4} {4, 1}
    // {4, 3} {3, 4}
    // {3, 2} {2, 3}
    // {2, 4} {4, 2}
    //
    // Edges after edge swap:
    // {0, 2} -> {0, 4}
    // {2, 4} -> {2, 2}

    DHBGraph G(this->Ghouse);
    ASSERT_EQ(0u, G.numberOfSelfLoops());

    NetworKit::node source_a = 0;
    NetworKit::node target_a = 2;
    NetworKit::node source_b = 2;
    NetworKit::node target_b = 4;

    ASSERT_TRUE(G.hasEdge(source_a, target_a));
    ASSERT_TRUE(G.hasEdge(source_b, target_b));

    ASSERT_FALSE(G.hasEdge(source_a, target_b));
    ASSERT_FALSE(G.hasEdge(source_b, target_a));

    float const weight_a = 5.f;
    float const weight_b = 3.f;

    if (G.isWeighted()) {
        G.setWeight(source_a, target_a, weight_a);
        G.setWeight(source_b, target_b, weight_b);
    }

    G.indexEdges();

    edgeid const id_origin_a = G.edgeId(source_a, target_a);
    edgeid const id_origin_b = G.edgeId(source_b, target_b);

    G.swapEdge(source_a, target_a, source_b, target_b);

    ASSERT_TRUE(G.hasEdge(source_a, target_b));
    ASSERT_TRUE(G.hasEdge(source_b, target_a));

    if (G.isWeighted()) {
        ASSERT_FLOAT_EQ(5.f, G.weight(source_a, target_b));
        ASSERT_FLOAT_EQ(3.f, G.weight(source_b, target_a));
    }

    ASSERT_EQ(id_origin_a, G.edgeId(source_a, target_b));
    ASSERT_EQ(id_origin_b, G.edgeId(source_b, target_a));

    ASSERT_EQ(1u, G.numberOfSelfLoops());
}
