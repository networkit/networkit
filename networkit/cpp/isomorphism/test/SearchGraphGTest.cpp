/*
 * SearchGraphGTest.cpp
 *
 *  Created on: Aug 12, 2026
 *      Author: Alexandra
 */

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/EdgeListReader.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"
#include "../SearchGraph.hpp"

namespace NetworKit {

using IsomorphismDetails::SearchGraph;

namespace {

/// The distinct out-neighbours, or with @a in the in-neighbours, of @a u other than @a u, in
/// ascending order.
std::vector<node> simpleNeighbors(const Graph &G, node u, bool in) {
    std::vector<node> neighbors;
    const auto collect = [&](node v) {
        if (v != u)
            neighbors.push_back(v);
    };
    if (in)
        G.forInNeighborsOf(u, collect);
    else
        G.forNeighborsOf(u, collect);
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    return neighbors;
}

void expectSlicesAreSimpleNeighborhoods(const Graph &G, const SearchGraph &SG) {
    G.forNodes([&](node u) {
        const std::vector<node> expectedOut = simpleNeighbors(G, u, /* in = */ false);
        const std::vector<node> expectedIn = simpleNeighbors(G, u, /* in = */ true);

        EXPECT_EQ(std::vector<node>(SG.outBegin(u), SG.outEnd(u)), expectedOut);
        EXPECT_EQ(std::vector<node>(SG.inBegin(u), SG.inEnd(u)), expectedIn);
        EXPECT_EQ(SG.outDegree(u), expectedOut.size());
        EXPECT_EQ(SG.inDegree(u), expectedIn.size());
    });
}

/// Every arc must keep the label of its own edge. A collapsed run of parallel arcs keeps one of the
/// given labels, so the check tests membership.
void expectLabelsPairWithHeads(const Graph &G, const std::vector<index> &edgeLabels,
                               const SearchGraph &SG) {
    const auto given = IsomorphismTest::edgeLabelsByPair(G, edgeLabels);

    G.forNodes([&](node u) {
        const index *outLabels = SG.outLabelBegin(u);
        ASSERT_NE(outLabels, nullptr);

        index offset = 0;
        for (const node *v = SG.outBegin(u); v != SG.outEnd(u); ++v, ++offset) {
            const std::vector<index> &wanted = IsomorphismTest::labelsOfPair(given, u, *v);
            EXPECT_NE(std::find(wanted.begin(), wanted.end(), outLabels[offset]), wanted.end())
                << "out-arc " << u << " -> " << *v << " carries a label nobody gave it";
            EXPECT_EQ(SG.edgeLabel(u, *v), outLabels[offset])
                << "edgeLabel() and the out-slice disagree about " << u << " -> " << *v;
        }

        const index *inLabels = SG.inLabelBegin(u);
        ASSERT_NE(inLabels, nullptr);

        offset = 0;
        for (const node *v = SG.inBegin(u); v != SG.inEnd(u); ++v, ++offset) {
            const std::vector<index> &wanted = IsomorphismTest::labelsOfPair(given, *v, u);
            EXPECT_NE(std::find(wanted.begin(), wanted.end(), inLabels[offset]), wanted.end())
                << "in-arc " << *v << " -> " << u << " carries a label nobody gave it";
        }
    });
}

} // namespace

class SearchGraphGTest : public testing::Test {};

TEST_F(SearchGraphGTest, testSearchGraphCSR) {

    Graph G = Graph(10);

    G.addEdge(3, 5);
    G.addEdge(9, 1);
    G.addEdge(5, 7);
    G.addEdge(1, 5);
    G.addEdge(5, 9);
    G.addEdge(5, 5);
    G.addEdge(2, 4);

    G.removeNode(0);
    G.removeNode(6);
    G.removeEdge(2, 4);

    SearchGraph SG(G, false);

    expectSlicesAreSimpleNeighborhoods(G, SG);

    // Even nodes never have incident edges, so none should be found
    G.forNodes([&](node u) {
        EXPECT_EQ(SG.outDegree(u), SG.inDegree(u));
        if (u % 2 == 0) {
            for (node v = 0; v < G.upperNodeIdBound(); v++) {
                EXPECT_FALSE(SG.hasEdge(u, v));
                EXPECT_FALSE(SG.hasEdge(v, u));
            }
        }
    });

    // The snapshot drops the self-loop at node 5.
    G.forEdges([&](node u, node v) {
        if (u == v) {
            EXPECT_FALSE(SG.hasEdge(u, v));
        } else {
            EXPECT_TRUE(SG.hasEdge(u, v));
            EXPECT_TRUE(SG.hasEdge(v, u));
        }
    });

    EXPECT_EQ(std::vector<node>(SG.inBegin(5), SG.inEnd(5)), (std::vector<node>{1, 3, 7, 9}));
    EXPECT_EQ(std::vector<node>(SG.outBegin(5), SG.outEnd(5)), (std::vector<node>{1, 3, 7, 9}));

    // G still has the self-loop, so the snapshot degree is one lower than the Graph degree
    EXPECT_EQ(G.degreeOut(5), 5);
    EXPECT_EQ(SG.outDegree(5), 4);
    EXPECT_FALSE(SG.isDirected());

    Graph D = Graph(8, false, true);

    D.addEdge(0, 5);
    D.addEdge(0, 4);
    D.addEdge(0, 3);
    D.addEdge(0, 2);
    D.addEdge(0, 1);
    D.addEdge(0, 0);

    D.removeNode(4);
    D.removeNode(7);
    D.removeEdge(0, 5);

    SearchGraph SD(D, false);

    expectSlicesAreSimpleNeighborhoods(D, SD);

    EXPECT_TRUE(SD.isDirected());
    EXPECT_NE(SD.inDegree(0), SD.outDegree(0));

    D.forEdges([&](node u, node v) {
        if (u == v) {
            EXPECT_FALSE(SD.hasEdge(u, v));
        } else {
            EXPECT_TRUE(SD.hasEdge(u, v));
            EXPECT_FALSE(SD.hasEdge(v, u));
        }
    });

    EXPECT_FALSE(SD.hasEdge(1, 2));
    EXPECT_FALSE(SD.hasEdge(0, 4));
    EXPECT_FALSE(SD.hasEdge(0, 5));

    EXPECT_EQ(std::vector<node>(SD.inBegin(0), SD.inEnd(0)), (std::vector<node>{}));
    EXPECT_EQ(std::vector<node>(SD.outBegin(0), SD.outEnd(0)), (std::vector<node>{1, 2, 3}));
    EXPECT_EQ(std::vector<node>(SD.inBegin(1), SD.inEnd(1)), (std::vector<node>{0}));
    EXPECT_EQ(std::vector<node>(SD.outBegin(1), SD.outEnd(1)), (std::vector<node>{}));

    // Deleted or isolated nodes have empty slice
    EXPECT_EQ(SD.inBegin(4), SD.inEnd(4));
    EXPECT_EQ(SD.inBegin(6), SD.inEnd(6));

    EXPECT_EQ(SD.numberOfNodes(), 6);
    EXPECT_EQ(SD.upperNodeIdBound(), 8);
}

TEST_F(SearchGraphGTest, testSearchGraphAdjMatrix) {

    // Graph G is directed and needs two 64-bit words per row
    Graph G = Graph(75, false, true);

    G.addEdge(1, 20);
    G.addEdge(42, 7);
    G.addEdge(50, 74);
    G.addEdge(38, 19);
    G.addEdge(25, 62);

    G.removeNode(38);
    G.removeEdge(25, 62);

    SearchGraph SG(G, true);

    G.forEdges([&](node u, node v) {
        EXPECT_TRUE(SG.hasEdge(u, v));
        EXPECT_FALSE(SG.hasEdge(v, u));
    });

    G.forNodes([&](node u) {
        if ((u != 1) && (u != 7) && (u != 20) && (u != 42) && (u != 50) && (u != 74)) {
            for (node v = 0; v < G.upperNodeIdBound(); v++) {
                EXPECT_FALSE(SG.hasEdge(u, v));
                EXPECT_FALSE(SG.hasEdge(v, u));
            }
        }
    });

    // Graph H is directed but has edge (v,u) for every edge (u,v) and needs one 64-bit word per row
    EdgeListReader reader('\t', 0, "#", true, true);
    Graph H = reader.read("input/example.edgelist");

    SearchGraph SH(H, true);

    H.forEdges([&](node u, node v) {
        EXPECT_TRUE(SH.hasEdge(u, v));
        EXPECT_TRUE(SH.hasEdge(v, u));
    });
}

TEST_F(SearchGraphGTest, testMultiEdgesAndLoopsCollapsed) {

    // Both hasEdge() backends must agree edge for edge. The id bound exceeds 64, so a matrix row
    // spans two words.
    for (bool directed : {false, true}) {
        SCOPED_TRACE(testing::Message() << "directed=" << directed);
        Graph G = Graph(70, false, directed);

        G.addEdge(1, 2);
        G.addEdge(1, 2); // parallel
        G.addEdge(2, 1); // the reverse; under `directed` a genuinely different edge
        G.addEdge(3, 4);
        G.addEdge(5, 5);
        G.addEdge(5, 5); // a repeated self-loop collapses to nothing at all
        G.addEdge(64, 65);
        G.addEdge(64, 65); // second word of each matrix row
        G.addEdge(65, 64);
        G.addEdge(69, 69);
        G.addEdge(7, 66); // across the word boundary
        G.addEdge(66, 7);

        if (directed) {
            ASSERT_EQ(G.degreeOut(1), 2);
            ASSERT_EQ(G.degreeIn(2), 2);
        } else {
            ASSERT_EQ(G.degreeOut(1), 3);
        }
        ASSERT_EQ(G.degreeOut(5), 2);

        SearchGraph S_CSR(G, false);
        SearchGraph S_Adj(G, true);

        expectSlicesAreSimpleNeighborhoods(G, S_CSR);
        expectSlicesAreSimpleNeighborhoods(G, S_Adj);

        for (const SearchGraph *SG : {&S_CSR, &S_Adj}) {
            EXPECT_EQ(std::vector<node>(SG->outBegin(1), SG->outEnd(1)), (std::vector<node>{2}));
            EXPECT_EQ(std::vector<node>(SG->outBegin(2), SG->outEnd(2)), (std::vector<node>{1}));
            EXPECT_EQ(SG->outDegree(1), 1);
            EXPECT_EQ(SG->outDegree(2), 1);

            EXPECT_EQ(std::vector<node>(SG->inBegin(2), SG->inEnd(2)), (std::vector<node>{1}));
            EXPECT_EQ(std::vector<node>(SG->inBegin(1), SG->inEnd(1)), (std::vector<node>{2}));

            EXPECT_TRUE(SG->hasEdge(1, 2));
            EXPECT_TRUE(SG->hasEdge(2, 1));
            EXPECT_TRUE(SG->hasEdge(3, 4));

            // Collapsing must not symmetrize anything.
            EXPECT_EQ(SG->hasEdge(4, 3), !directed);
            EXPECT_FALSE(SG->hasEdge(1, 3));

            EXPECT_EQ(SG->outDegree(5), 0);
            EXPECT_EQ(SG->inDegree(5), 0);
            EXPECT_FALSE(SG->hasEdge(5, 5));
            EXPECT_FALSE(SG->hasEdge(69, 69));
        }

        for (node u = 0; u < G.upperNodeIdBound(); u++) {
            for (node v = 0; v < G.upperNodeIdBound(); v++) {
                EXPECT_EQ(S_CSR.hasEdge(u, v), S_Adj.hasEdge(u, v)) << "u=" << u << " v=" << v;
            }
            EXPECT_FALSE(S_CSR.hasEdge(u, u));
            EXPECT_FALSE(S_Adj.hasEdge(u, u));
        }
    }
}

TEST_F(SearchGraphGTest, testMatrixFallsBackForLargeIdBound) {

    // upperNodeIdBound() sets the matrix size, and declining the matrix must not change hasEdge().
    Graph big = Graph(30000);
    big.addEdge(1, 2);

    SearchGraph SBig(big, true);

    EXPECT_FALSE(SBig.hasAdjacencyMatrix());
    EXPECT_TRUE(SBig.hasEdge(1, 2));
    EXPECT_TRUE(SBig.hasEdge(2, 1));
    EXPECT_FALSE(SBig.hasEdge(1, 3));
    EXPECT_FALSE(SBig.hasEdge(29998, 29999));

    // The matrix may take 64 MiB. With rows of whole 64-bit words, that limit falls at 23168 ids.
    EXPECT_TRUE(SearchGraph(Graph(23168), true).hasAdjacencyMatrix());
    EXPECT_FALSE(SearchGraph(Graph(23169), true).hasAdjacencyMatrix());

    // removeNode() never lowers the id bound.
    Graph carved = Graph(30000);
    carved.addEdge(0, 1);
    carved.addEdge(1, 2);
    for (node u = 3; u < 30000; ++u)
        carved.removeNode(u);

    ASSERT_EQ(carved.numberOfNodes(), 3);
    ASSERT_EQ(carved.upperNodeIdBound(), 30000);

    SearchGraph SCarved(carved, true);

    EXPECT_FALSE(SCarved.hasAdjacencyMatrix());
    EXPECT_TRUE(SCarved.hasEdge(0, 1));
    EXPECT_TRUE(SCarved.hasEdge(1, 2));
    EXPECT_FALSE(SCarved.hasEdge(0, 2));

    // The warning recommends compacting the ids.
    Graph compacted =
        GraphTools::getCompactedGraph(carved, GraphTools::getContinuousNodeIds(carved));

    ASSERT_EQ(compacted.upperNodeIdBound(), 3);

    SearchGraph SCompacted(compacted, true);

    EXPECT_TRUE(SCompacted.hasAdjacencyMatrix());
    expectSlicesAreSimpleNeighborhoods(compacted, SCompacted);

    Graph ordinary = Graph(75);
    ordinary.addEdge(1, 20);
    ordinary.removeNode(38);

    SearchGraph SOrdinary(ordinary, true);

    EXPECT_TRUE(SOrdinary.hasAdjacencyMatrix());
    EXPECT_TRUE(SOrdinary.hasEdge(1, 20));

    EXPECT_FALSE(SearchGraph(ordinary, false).hasAdjacencyMatrix());
}

TEST_F(SearchGraphGTest, testHasNode) {

    Graph G = Graph(10);
    G.addEdge(0, 1);
    G.removeNode(4);
    G.removeNode(9); // the last id: the bound must not shrink with it

    SearchGraph SG(G, true);

    ASSERT_EQ(SG.upperNodeIdBound(), 10);
    ASSERT_EQ(SG.numberOfNodes(), 8);

    for (node u = 0; u < SG.upperNodeIdBound(); ++u) {
        EXPECT_EQ(SG.hasNode(u), G.hasNode(u)) << "disagreement about node " << u;
    }

    // Node 7 is isolated and node 4 is removed, so both have empty slices.
    EXPECT_EQ(SG.outDegree(7), SG.outDegree(4));
    EXPECT_TRUE(SG.hasNode(7));
    EXPECT_FALSE(SG.hasNode(4));

    SearchGraph SEmpty(Graph(0), true);
    EXPECT_EQ(SEmpty.numberOfNodes(), 0);
    EXPECT_EQ(SEmpty.upperNodeIdBound(), 0);

    Graph allRemoved = Graph(5);
    for (node u = 0; u < 5; ++u)
        allRemoved.removeNode(u);

    SearchGraph SAllRemoved(allRemoved, true);

    ASSERT_EQ(SAllRemoved.numberOfNodes(), 0);
    ASSERT_EQ(SAllRemoved.upperNodeIdBound(), 5);
    for (node u = 0; u < 5; ++u) {
        EXPECT_FALSE(SAllRemoved.hasNode(u));
    }

    SearchGraph SFull(Graph(3), true);
    for (node u = 0; u < 3; ++u) {
        EXPECT_TRUE(SFull.hasNode(u));
    }
}

TEST_F(SearchGraphGTest, testMaxDegree) {

    // Parallel edges and loops make the maximum differ from GraphTools::maxDegree().
    Graph G = Graph(5);
    G.addEdge(0, 1);
    G.addEdge(0, 1); // parallel
    G.addEdge(0, 1); // parallel
    G.addEdge(0, 2);
    G.addEdge(3, 3); // self-loop
    G.addEdge(3, 4);

    SearchGraph SG(G, true);

    ASSERT_EQ(G.degree(0), 4);
    EXPECT_EQ(SG.outDegree(0), 2);
    EXPECT_EQ(SG.maxOutDegree(), 2);
    EXPECT_EQ(SG.maxInDegree(), 2); // undirected: mirrors maxOutDegree()

    EXPECT_NE(SG.maxOutDegree(), GraphTools::maxDegree(G));

    // The two maxima come from different nodes.
    Graph D = Graph(4, false, true);
    D.addEdge(0, 1);
    D.addEdge(0, 2);
    D.addEdge(0, 3); // out-degree 3, in-degree 0
    D.addEdge(1, 3);
    D.addEdge(2, 3); // node 3 has in-degree 3, out-degree 0

    SearchGraph SD(D, true);

    EXPECT_EQ(SD.maxOutDegree(), 3);
    EXPECT_EQ(SD.maxInDegree(), 3);
    EXPECT_EQ(SD.outDegree(3), 0);
    EXPECT_EQ(SD.inDegree(3), 3);

    EXPECT_EQ(SearchGraph(Graph(0), false).maxOutDegree(), 0);
    EXPECT_EQ(SearchGraph(Graph(0), false).maxInDegree(), 0);
    EXPECT_EQ(SearchGraph(Graph(5), false).maxOutDegree(), 0);
}

TEST_F(SearchGraphGTest, testEdgeLabelsStayWithTheirArcs) {

    // Node 9 gets its neighbours in descending order, so the sort must move every label.
    for (bool directed : {false, true}) {
        SCOPED_TRACE(testing::Message() << "directed=" << directed);
        const IsomorphismTest::LabelledGraph labelled = IsomorphismTest::labelledGraphOf(
            10, {{9, 6, 60}, {9, 4, 40}, {9, 2, 20}, {9, 0, 5}, {8, 7, 70}, {8, 1, 10}, {3, 5, 50}},
            directed);

        ASSERT_EQ(std::vector<node>(labelled.G.neighborRange(9).begin(),
                                    labelled.G.neighborRange(9).end()),
                  (std::vector<node>{6, 4, 2, 0}));

        for (bool buildMatrix : {false, true}) {
            const SearchGraph SG(labelled.G, buildMatrix, labelled.edgeLabels);

            EXPECT_TRUE(SG.hasEdgeLabels());
            expectSlicesAreSimpleNeighborhoods(labelled.G, SG);
            expectLabelsPairWithHeads(labelled.G, labelled.edgeLabels, SG);

            EXPECT_EQ(std::vector<node>(SG.outBegin(9), SG.outEnd(9)),
                      (std::vector<node>{0, 2, 4, 6}));
            EXPECT_EQ(std::vector<index>(SG.outLabelBegin(9), SG.outLabelBegin(9) + 4),
                      (std::vector<index>{5, 20, 40, 60}));

            // The reverse of a directed arc has no label.
            EXPECT_EQ(SG.edgeLabel(9, 7), none);
            EXPECT_EQ(SG.edgeLabel(6, 9), directed ? none : index{60});

            EXPECT_FALSE(SG.collapsedLabelledEdges());
        }
    }

    // The compaction must move the labels too. Parallel arcs agree on their labels here.
    for (bool directed : {false, true}) {
        SCOPED_TRACE(testing::Message() << "directed=" << directed);
        const IsomorphismTest::LabelledGraph labelled =
            IsomorphismTest::labelledGraphOf(6,
                                             {{4, 3, 30},
                                              {4, 3, 30},
                                              {4, 1, 10},
                                              {4, 1, 10},
                                              {4, 1, 10},
                                              {2, 2, 99},
                                              {2, 2, 98},
                                              {5, 0, 50}},
                                             directed);

        const SearchGraph SG(labelled.G, /* buildMatrix = */ true, labelled.edgeLabels);

        expectSlicesAreSimpleNeighborhoods(labelled.G, SG);
        expectLabelsPairWithHeads(labelled.G, labelled.edgeLabels, SG);

        EXPECT_EQ(std::vector<node>(SG.outBegin(4), SG.outEnd(4)), (std::vector<node>{1, 3}));
        EXPECT_EQ(SG.edgeLabel(4, 1), 10u);
        EXPECT_EQ(SG.edgeLabel(4, 3), 30u);

        EXPECT_EQ(SG.outDegree(2), 0u);
        EXPECT_EQ(SG.edgeLabel(2, 2), none);
        EXPECT_FALSE(SG.collapsedLabelledEdges());
    }

    // The two arcs of a directed mutual pair carry independent labels.
    const IsomorphismTest::LabelledGraph mutual = IsomorphismTest::labelledGraphOf(
        3, {{0, 1, 7}, {1, 0, 8}, {1, 2, 9}}, /* directed = */ true);

    const SearchGraph SM(mutual.G, /* buildMatrix = */ true, mutual.edgeLabels);

    EXPECT_EQ(SM.edgeLabel(0, 1), 7u);
    EXPECT_EQ(SM.edgeLabel(1, 0), 8u);
    EXPECT_EQ(SM.edgeLabel(1, 2), 9u);
    EXPECT_EQ(SM.edgeLabel(2, 1), none);
    EXPECT_FALSE(SM.collapsedLabelledEdges());

    expectLabelsPairWithHeads(mutual.G, mutual.edgeLabels, SM);

    // The in-slices must carry the same labels, arc for arc.
    ASSERT_EQ(SM.inDegree(1), 1u);
    EXPECT_EQ(*SM.inBegin(1), 0u);
    EXPECT_EQ(*SM.inLabelBegin(1), 7u);

    ASSERT_EQ(SM.inDegree(0), 1u);
    EXPECT_EQ(*SM.inBegin(0), 1u);
    EXPECT_EQ(*SM.inLabelBegin(0), 8u);
}

TEST_F(SearchGraphGTest, testCollapsedLabelledEdges) {

    for (bool directed : {false, true}) {
        SCOPED_TRACE(testing::Message() << "directed=" << directed);
        const IsomorphismTest::LabelledGraph plain =
            IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}}, directed);
        const IsomorphismTest::LabelledGraph agreeing = IsomorphismTest::labelledGraphOf(
            4, {{0, 1, 1}, {0, 1, 1}, {1, 2, 2}, {2, 3, 1}}, directed);
        const IsomorphismTest::LabelledGraph disagreeing = IsomorphismTest::labelledGraphOf(
            4, {{0, 1, 1}, {0, 1, 4}, {1, 2, 2}, {2, 3, 1}}, directed);

        EXPECT_FALSE(SearchGraph(plain.G, false, plain.edgeLabels).collapsedLabelledEdges())
            << "nothing was collapsed at all";

        EXPECT_FALSE(SearchGraph(agreeing.G, false, agreeing.edgeLabels).collapsedLabelledEdges())
            << "collapsing equal labels is lossless";

        EXPECT_TRUE(
            SearchGraph(disagreeing.G, false, disagreeing.edgeLabels).collapsedLabelledEdges())
            << "a label was thrown away";

        const SearchGraph unlabelled(disagreeing.G, false);
        EXPECT_FALSE(unlabelled.hasEdgeLabels());
        EXPECT_FALSE(unlabelled.collapsedLabelledEdges());
        EXPECT_EQ(unlabelled.edgeLabel(0, 1), none);
        EXPECT_EQ(unlabelled.outLabelBegin(0), nullptr);
    }
}

} // namespace NetworKit
