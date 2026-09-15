/*
 * SearchGraphGTest.cpp
 *
 *  Created on: Aug 12, 2026
 *      Author: Alexandra
 */

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"
#include "../SearchGraph.hpp"

namespace NetworKit {

namespace {

/// What SearchGraph promises an out-slice contains: the distinct out-neighbours of @a u, without
/// @a u itself, sorted ascending. Not the same as G.degreeOut(u) once G has loops or multi-edges.
std::vector<node> simpleOutNeighbors(const Graph &G, node u) {
    std::vector<node> neighbors;
    G.forNeighborsOf(u, [&](node v) {
        if (v != u)
            neighbors.push_back(v);
    });
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    return neighbors;
}

std::vector<node> simpleInNeighbors(const Graph &G, node u) {
    std::vector<node> neighbors;
    G.forInNeighborsOf(u, [&](node v) {
        if (v != u)
            neighbors.push_back(v);
    });
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    return neighbors;
}

/// Every slice must match the simple-graph neighbourhood exactly, in both directions.
void expectSlicesAreSimpleNeighborhoods(const Graph &G, const IsomorphismDetails::SearchGraph &SG) {
    G.forNodes([&](node u) {
        const std::vector<node> expectedOut = simpleOutNeighbors(G, u);
        const std::vector<node> expectedIn = simpleInNeighbors(G, u);

        EXPECT_EQ(std::vector<node>(SG.outBegin(u), SG.outEnd(u)), expectedOut);
        EXPECT_EQ(std::vector<node>(SG.inBegin(u), SG.inEnd(u)), expectedIn);
        EXPECT_EQ(SG.outDegree(u), expectedOut.size());
        EXPECT_EQ(SG.inDegree(u), expectedIn.size());
    });
}

/// Every arc must still be holding the label its own edge was given.
///
/// This is the invariant the whole edge-label layout rests on: a label is identified by the offset
/// it sits at, so the sort and the compaction have to move the label array in lockstep with the
/// head array. Getting that wrong breaks nothing loudly - every arc simply ends up with some other
/// arc's label - which is why it is checked directly rather than through anything that uses it.
///
/// Both ways of reading a label are checked against each other and against the input: the
/// @ref SearchGraph::outLabelBegin() slice, which the search will walk, and
/// @ref SearchGraph::edgeLabel(), which binary-searches for one arc. Membership rather than
/// equality against the input, because a collapsed run of parallel arcs legitimately keeps only one
/// of the labels the caller gave that pair; without parallel edges there is only one to keep.
void expectLabelsPairWithHeads(const Graph &G, const std::vector<index> &edgeLabels,
                               const IsomorphismDetails::SearchGraph &SG) {
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

    // The CSR backend answers the same questions in both directednesses, so both are one test.
    // Undirected first: every edge is reachable from either end.
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

    IsomorphismDetails::SearchGraph SG = IsomorphismDetails::SearchGraph(G, false);

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

    // All present edges must be found (oriented both ways). The self-loop at node 5 is dropped,
    // so it must NOT be found - both semantics only ever constrain pairs of distinct nodes.
    G.forEdges([&](node u, node v) {
        if (u == v) {
            EXPECT_FALSE(SG.hasEdge(u, v));
        } else {
            EXPECT_TRUE(SG.hasEdge(u, v));
            EXPECT_TRUE(SG.hasEdge(v, u));
        }
    });

    // Check in and out slice for node 5. The self-loop is not part of either.
    EXPECT_EQ(std::vector<node>(SG.inBegin(5), SG.inEnd(5)), (std::vector<node>{1, 3, 7, 9}));
    EXPECT_EQ(std::vector<node>(SG.outBegin(5), SG.outEnd(5)), (std::vector<node>{1, 3, 7, 9}));
    EXPECT_EQ(std::vector<node>(SG.inBegin(5), SG.inEnd(5)),
              std::vector<node>(SG.outBegin(5), SG.outEnd(5)));

    // G still has the self-loop, so the snapshot degree is one lower than the Graph degree
    EXPECT_EQ(G.degreeOut(5), 5);
    EXPECT_EQ(SG.outDegree(5), 4);
    EXPECT_FALSE(SG.isDirected());

    // Directed now. The same rules hold, except that an arc is reachable from one end only. So a
    // node's two degrees come apart and hasEdge() stops being symmetric.
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

    IsomorphismDetails::SearchGraph SD = IsomorphismDetails::SearchGraph(D, false);

    expectSlicesAreSimpleNeighborhoods(D, SD);

    EXPECT_TRUE(SD.isDirected());
    EXPECT_NE(SD.inDegree(0), SD.outDegree(0));

    // All present edges must be found, except the self loop at node 0, which is dropped
    D.forEdges([&](node u, node v) {
        if (u == v) {
            EXPECT_FALSE(SD.hasEdge(u, v));
        } else {
            EXPECT_TRUE(SD.hasEdge(u, v));
            EXPECT_FALSE(SD.hasEdge(v, u));
        }
    });

    // Non-present edges must not be found
    EXPECT_FALSE(SD.hasEdge(1, 2));
    EXPECT_FALSE(SD.hasEdge(0, 4));
    EXPECT_FALSE(SD.hasEdge(0, 5));

    // Check in and out slices. Node 0's self-loop appears in neither.
    EXPECT_EQ(std::vector<node>(SD.inBegin(0), SD.inEnd(0)), (std::vector<node>{}));
    EXPECT_EQ(std::vector<node>(SD.outBegin(0), SD.outEnd(0)), (std::vector<node>{1, 2, 3}));
    EXPECT_EQ(std::vector<node>(SD.inBegin(1), SD.inEnd(1)), (std::vector<node>{0}));
    EXPECT_EQ(std::vector<node>(SD.outBegin(1), SD.outEnd(1)), (std::vector<node>{}));

    // Deleted or isolated nodes have empty slice
    EXPECT_EQ(SD.inBegin(4), SD.inEnd(4));
    EXPECT_EQ(SD.inBegin(6), SD.inEnd(6));

    // Both snapshots keep the node counts the Graph had, removed ids included in the bound.
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

    IsomorphismDetails::SearchGraph SG = IsomorphismDetails::SearchGraph(G, true);

    // All present edges must be found
    G.forEdges([&](node u, node v) {
        EXPECT_TRUE(SG.hasEdge(u, v));
        EXPECT_FALSE(SG.hasEdge(v, u));
    });

    // Non-present edges must not be found
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

    IsomorphismDetails::SearchGraph SH = IsomorphismDetails::SearchGraph(H, true);

    H.forEdges([&](node u, node v) {
        EXPECT_TRUE(SH.hasEdge(u, v));
        EXPECT_TRUE(SH.hasEdge(v, u));
    });
}

TEST_F(SearchGraphGTest, testMultiEdgesAndLoopsCollapsed) {

    // Graph::addEdge() permits parallel edges and self-loops by default, so a snapshot has to
    // collapse them: a repeated neighbour would make a search enumerate the same candidate twice
    // and would inflate the degrees every feasibility rule prunes on, while a self-loop is never
    // usable because both semantics only constrain pairs of *distinct* nodes.
    //
    // Both directednesses and both hasEdge() backends have to agree about all of that, which is
    // why they are one test: the CSR drops self-loops, so the matrix has to leave the diagonal
    // clear too. This is also the graph the two backends are compared on edge for edge, because
    // an ordinary graph has neither loops nor parallel edges to catch a disagreement with. The id
    // bound reaches past 64 so the matrix spans more than one word per row.
    for (bool directed : {false, true}) {
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

        IsomorphismDetails::SearchGraph S_CSR = IsomorphismDetails::SearchGraph(G, false);
        IsomorphismDetails::SearchGraph S_Adj = IsomorphismDetails::SearchGraph(G, true);

        expectSlicesAreSimpleNeighborhoods(G, S_CSR);
        expectSlicesAreSimpleNeighborhoods(G, S_Adj);

        for (const IsomorphismDetails::SearchGraph *SG : {&S_CSR, &S_Adj}) {
            EXPECT_EQ(std::vector<node>(SG->outBegin(1), SG->outEnd(1)), (std::vector<node>{2}))
                << "directed=" << directed;
            EXPECT_EQ(std::vector<node>(SG->outBegin(2), SG->outEnd(2)), (std::vector<node>{1}))
                << "directed=" << directed;
            EXPECT_EQ(SG->outDegree(1), 1) << "directed=" << directed;
            EXPECT_EQ(SG->outDegree(2), 1) << "directed=" << directed;

            EXPECT_EQ(std::vector<node>(SG->inBegin(2), SG->inEnd(2)), (std::vector<node>{1}))
                << "directed=" << directed;
            EXPECT_EQ(std::vector<node>(SG->inBegin(1), SG->inEnd(1)), (std::vector<node>{2}))
                << "directed=" << directed;

            EXPECT_TRUE(SG->hasEdge(1, 2)) << "directed=" << directed;
            EXPECT_TRUE(SG->hasEdge(2, 1)) << "directed=" << directed;
            EXPECT_TRUE(SG->hasEdge(3, 4)) << "directed=" << directed;

            // Collapsing must not symmetrize anything: 3-4 was only ever added one way round.
            EXPECT_EQ(SG->hasEdge(4, 3), !directed) << "directed=" << directed;
            EXPECT_FALSE(SG->hasEdge(1, 3)) << "directed=" << directed;

            // Self-loops leave no trace, in either the slices or the degrees.
            EXPECT_EQ(SG->outDegree(5), 0) << "directed=" << directed;
            EXPECT_EQ(SG->inDegree(5), 0) << "directed=" << directed;
            EXPECT_FALSE(SG->hasEdge(5, 5)) << "directed=" << directed;
            EXPECT_FALSE(SG->hasEdge(69, 69)) << "directed=" << directed;
        }

        // The two backends must not drift apart anywhere, on or off the diagonal.
        for (node u = 0; u < G.upperNodeIdBound(); u++) {
            for (node v = 0; v < G.upperNodeIdBound(); v++) {
                EXPECT_EQ(S_CSR.hasEdge(u, v), S_Adj.hasEdge(u, v))
                    << "directed=" << directed << " u=" << u << " v=" << v;
            }
            // No node is ever its own neighbour, whichever backend answers
            EXPECT_FALSE(S_CSR.hasEdge(u, u)) << "directed=" << directed;
            EXPECT_FALSE(S_Adj.hasEdge(u, u)) << "directed=" << directed;
        }
    }
}

TEST_F(SearchGraphGTest, testSlicesAreStrictlyAscending) {

    // Strictly ascending proves sortedness and dedup at once, and hasEdge() binary-searches the
    // slice, so this is the invariant the whole class rests on.
    METISGraphReader reader;
    Graph G = reader.read("input/karate.graph");

    IsomorphismDetails::SearchGraph SG = IsomorphismDetails::SearchGraph(G, false);

    G.forNodes([&](node u) {
        for (const node *it = SG.outBegin(u); it != SG.outEnd(u); ++it) {
            EXPECT_NE(*it, u) << "node " << u << " is its own neighbour";
            if (it + 1 != SG.outEnd(u)) {
                EXPECT_LT(*it, *(it + 1)) << "out-slice of " << u << " is not strictly ascending";
            }
        }
        for (const node *it = SG.inBegin(u); it != SG.inEnd(u); ++it) {
            EXPECT_NE(*it, u);
            if (it + 1 != SG.inEnd(u)) {
                EXPECT_LT(*it, *(it + 1)) << "in-slice of " << u << " is not strictly ascending";
            }
        }
    });
}

TEST_F(SearchGraphGTest, testMatrixFallsBackForLargeIdBound) {

    // The matrix needs one bit per ordered pair of node *ids*, so it is sized by
    // upperNodeIdBound() rather than by how many nodes exist. Asking for it on anything but a
    // small pattern is a request the snapshot declines - and the decline has to be invisible,
    // because hasEdge() answers from the CSR either way.
    Graph big = Graph(30000);
    big.addEdge(1, 2);

    IsomorphismDetails::SearchGraph SBig = IsomorphismDetails::SearchGraph(big, true);

    EXPECT_FALSE(SBig.hasAdjacencyMatrix());
    EXPECT_TRUE(SBig.hasEdge(1, 2));
    EXPECT_TRUE(SBig.hasEdge(2, 1));
    EXPECT_FALSE(SBig.hasEdge(1, 3));
    EXPECT_FALSE(SBig.hasEdge(29998, 29999));

    // The trap this really guards: removeNode() clears a slot but never lowers the id bound, so
    // a three-node pattern carved out of a large graph still asks for the large graph's matrix.
    Graph carved = Graph(30000);
    carved.addEdge(0, 1);
    carved.addEdge(1, 2);
    for (node u = 3; u < 30000; ++u)
        carved.removeNode(u);

    ASSERT_EQ(carved.numberOfNodes(), 3);
    ASSERT_EQ(carved.upperNodeIdBound(), 30000);

    IsomorphismDetails::SearchGraph SCarved = IsomorphismDetails::SearchGraph(carved, true);

    EXPECT_FALSE(SCarved.hasAdjacencyMatrix());
    EXPECT_TRUE(SCarved.hasEdge(0, 1));
    EXPECT_TRUE(SCarved.hasEdge(1, 2));
    EXPECT_FALSE(SCarved.hasEdge(0, 2));

    // Compacting the ids is what the warning tells the caller to do, so it must actually work
    Graph compacted =
        GraphTools::getCompactedGraph(carved, GraphTools::getContinuousNodeIds(carved));

    ASSERT_EQ(compacted.upperNodeIdBound(), 3);

    IsomorphismDetails::SearchGraph SCompacted = IsomorphismDetails::SearchGraph(compacted, true);

    EXPECT_TRUE(SCompacted.hasAdjacencyMatrix());
    expectSlicesAreSimpleNeighborhoods(compacted, SCompacted);

    // A few removed nodes must not trip the guard - that is the ordinary case
    Graph ordinary = Graph(75);
    ordinary.addEdge(1, 20);
    ordinary.removeNode(38);

    IsomorphismDetails::SearchGraph SOrdinary = IsomorphismDetails::SearchGraph(ordinary, true);

    EXPECT_TRUE(SOrdinary.hasAdjacencyMatrix());
    EXPECT_TRUE(SOrdinary.hasEdge(1, 20));

    // Not requesting the matrix must leave it unbuilt whatever the size
    EXPECT_FALSE(IsomorphismDetails::SearchGraph(ordinary, false).hasAdjacencyMatrix());
}

TEST_F(SearchGraphGTest, testHasNode) {

    // The distinction hasNode() exists for: a removed id and an isolated node both have an empty
    // slice, so adjacency alone cannot tell them apart. A search enumerating "every unmapped
    // target node" that skips this check will map pattern nodes onto ids that are not nodes.
    Graph G = Graph(10);
    G.addEdge(0, 1);
    G.removeNode(4);
    G.removeNode(9); // the last id: the bound must not shrink with it

    IsomorphismDetails::SearchGraph SG = IsomorphismDetails::SearchGraph(G, true);

    ASSERT_EQ(SG.upperNodeIdBound(), 10);
    ASSERT_EQ(SG.numberOfNodes(), 8);

    for (node u = 0; u < SG.upperNodeIdBound(); ++u) {
        EXPECT_EQ(SG.hasNode(u), G.hasNode(u)) << "disagreement about node " << u;
    }

    // Node 7 exists and is isolated; node 4 does not exist. Identical slices, opposite answers.
    EXPECT_EQ(SG.outDegree(7), SG.outDegree(4));
    EXPECT_TRUE(SG.hasNode(7));
    EXPECT_FALSE(SG.hasNode(4));

    // The three degenerate graphs, where an `n == z` shortcut would read as a free optimisation
    // and be wrong: on the empty graph it makes hasNode() true for every id.
    IsomorphismDetails::SearchGraph SEmpty = IsomorphismDetails::SearchGraph(Graph(0), true);
    EXPECT_EQ(SEmpty.numberOfNodes(), 0);
    EXPECT_EQ(SEmpty.upperNodeIdBound(), 0);

    Graph allRemoved = Graph(5);
    for (node u = 0; u < 5; ++u)
        allRemoved.removeNode(u);

    IsomorphismDetails::SearchGraph SAllRemoved = IsomorphismDetails::SearchGraph(allRemoved, true);

    ASSERT_EQ(SAllRemoved.numberOfNodes(), 0);
    ASSERT_EQ(SAllRemoved.upperNodeIdBound(), 5);
    for (node u = 0; u < 5; ++u) {
        EXPECT_FALSE(SAllRemoved.hasNode(u));
    }

    // A graph where every id is a node is the other extreme, and must stay true throughout
    IsomorphismDetails::SearchGraph SFull = IsomorphismDetails::SearchGraph(Graph(3), true);
    for (node u = 0; u < 3; ++u) {
        EXPECT_TRUE(SFull.hasNode(u));
    }
}

TEST_F(SearchGraphGTest, testMaxDegree) {

    // The number callers compare a pattern degree against. Taken from the *compacted* slices, so
    // it disagrees with GraphTools::maxDegree() the moment there are parallel edges or a loop -
    // and it is this one that is right for a subgraph search.
    Graph G = Graph(5);
    G.addEdge(0, 1);
    G.addEdge(0, 1); // parallel
    G.addEdge(0, 1); // parallel
    G.addEdge(0, 2);
    G.addEdge(3, 3); // self-loop
    G.addEdge(3, 4);

    IsomorphismDetails::SearchGraph SG = IsomorphismDetails::SearchGraph(G, true);

    // Node 0 has four incident edges but only two distinct neighbours
    ASSERT_EQ(G.degree(0), 4);
    EXPECT_EQ(SG.outDegree(0), 2);
    EXPECT_EQ(SG.maxOutDegree(), 2);
    EXPECT_EQ(SG.maxInDegree(), 2); // undirected: mirrors maxOutDegree()

    count expected = 0;
    G.forNodes([&](node u) { expected = std::max(expected, SG.outDegree(u)); });
    EXPECT_EQ(SG.maxOutDegree(), expected);

    // The whole point: the naive number would have been 4, and pruning on it would reject valid
    // host nodes of degree 2 and 3 and silently lose real matches.
    EXPECT_NE(SG.maxOutDegree(), GraphTools::maxDegree(G));

    // Directed, where the two maxima come from different nodes and must not be conflated.
    Graph D = Graph(4, false, true);
    D.addEdge(0, 1);
    D.addEdge(0, 2);
    D.addEdge(0, 3); // out-degree 3, in-degree 0
    D.addEdge(1, 3);
    D.addEdge(2, 3); // node 3 has in-degree 3, out-degree 0

    IsomorphismDetails::SearchGraph SD = IsomorphismDetails::SearchGraph(D, true);

    EXPECT_EQ(SD.maxOutDegree(), 3);
    EXPECT_EQ(SD.maxInDegree(), 3);
    EXPECT_EQ(SD.outDegree(3), 0);
    EXPECT_EQ(SD.inDegree(3), 3);

    // No nodes and no edges both have to give 0 rather than reading past an empty array
    EXPECT_EQ(IsomorphismDetails::SearchGraph(Graph(0), false).maxOutDegree(), 0);
    EXPECT_EQ(IsomorphismDetails::SearchGraph(Graph(0), false).maxInDegree(), 0);
    EXPECT_EQ(IsomorphismDetails::SearchGraph(Graph(5), false).maxOutDegree(), 0);
}

TEST_F(SearchGraphGTest, testIntersectionSize) {

    // Checked against std::set_intersection rather than against hand-counted answers, so the test
    // does not encode the same off-by-one the implementation might.
    auto expectAgrees = [](std::vector<node> a, std::vector<node> b) {
        std::vector<node> common;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(common));
        EXPECT_EQ(IsomorphismDetails::intersectionSize(a.data(), a.data() + a.size(), b.data(),
                                                       b.data() + b.size()),
                  common.size());
    };

    expectAgrees({}, {});
    expectAgrees({1, 2, 3}, {});
    expectAgrees({}, {1, 2, 3});
    expectAgrees({1, 3, 5}, {2, 4, 6});    // disjoint
    expectAgrees({1, 2, 3}, {1, 2, 3});    // identical
    expectAgrees({1, 2, 3}, {3});          // last element only
    expectAgrees({1, 2, 3}, {1});          // first element only
    expectAgrees({0, 5, 9}, {5, 9, 11});   // partial overlap, unequal lengths
    expectAgrees({2}, {1, 2, 3, 4, 5, 6}); // one against many

    Aux::Random::setSeed(42, false);
    for (int trial = 0; trial < 50; ++trial) {
        std::vector<node> a, b;
        for (node u = 0; u < 40; ++u) {
            if (Aux::Random::probability() < 0.4)
                a.push_back(u);
            if (Aux::Random::probability() < 0.4)
                b.push_back(u);
        }
        expectAgrees(a, b);
    }
}

TEST_F(SearchGraphGTest, testCommonOutNeighbors) {

    // Every neighbourhood-cardinality pruning rule rests on this, so it must be the *distinct*
    // common neighbour count even when the input has parallel edges.
    Graph G = Graph(6);
    G.addEdge(0, 2);
    G.addEdge(0, 3);
    G.addEdge(0, 4);
    G.addEdge(1, 3);
    G.addEdge(1, 3); // parallel: must not double-count node 3
    G.addEdge(1, 4);
    G.addEdge(1, 5);

    IsomorphismDetails::SearchGraph SG = IsomorphismDetails::SearchGraph(G, false);

    // N(0) = {2,3,4}, N(1) = {3,4,5}, N(2) = {0}, N(3) = {0,1}, N(5) = {1}
    EXPECT_EQ(SG.commonOutNeighbors(0, 1), 2); // nodes 3 and 4, counted once despite the parallel
    EXPECT_EQ(SG.commonOutNeighbors(1, 0), 2); // symmetric when undirected
    EXPECT_EQ(SG.commonOutNeighbors(0, 5), 0); // {2,3,4} and {1} share nothing
    EXPECT_EQ(SG.commonOutNeighbors(2, 3), 1); // node 0
    EXPECT_EQ(SG.commonOutNeighbors(0, 0), SG.outDegree(0));
}

TEST_F(SearchGraphGTest, testEdgeLabelsStayWithTheirArcs) {

    // A label is identified by the offset it sits at, so every step that moves a head has to move
    // its label with it. Three steps can break that pairing, and all three are checked here,
    // because the symptom is the same in each case: an arc quietly ends up holding another arc's
    // label, and nothing else goes visibly wrong.

    // First the sort. Node 9 is given its neighbours in *descending* id order, so the scattered
    // order and the sorted order genuinely differ and std::sort has real work to do. Every edge
    // carries a different label, so a sort that moved heads without moving labels with them
    // cannot cancel out - each arc would visibly end up with another arc's label.
    for (bool directed : {false, true}) {
        const IsomorphismTest::LabelledGraph labelled = IsomorphismTest::labelledGraphOf(
            10, {{9, 6, 60}, {9, 4, 40}, {9, 2, 20}, {9, 0, 5}, {8, 7, 70}, {8, 1, 10}, {3, 5, 50}},
            directed);

        // The premise of the test: without this the sort is a no-op and nothing is being checked.
        ASSERT_EQ(std::vector<node>(labelled.G.neighborRange(9).begin(),
                                    labelled.G.neighborRange(9).end()),
                  (std::vector<node>{6, 4, 2, 0}))
            << "directed=" << directed;

        for (bool buildMatrix : {false, true}) {
            const IsomorphismDetails::SearchGraph SG(labelled.G, buildMatrix, labelled.edgeLabels);

            EXPECT_TRUE(SG.hasEdgeLabels()) << "directed=" << directed;
            expectSlicesAreSimpleNeighborhoods(labelled.G, SG);
            expectLabelsPairWithHeads(labelled.G, labelled.edgeLabels, SG);

            // Spelled out for one node, so a failure says what went wrong rather than only that
            // something did. The slice is sorted, the labels came in the opposite order.
            EXPECT_EQ(std::vector<node>(SG.outBegin(9), SG.outEnd(9)),
                      (std::vector<node>{0, 2, 4, 6}))
                << "directed=" << directed;
            EXPECT_EQ(std::vector<index>(SG.outLabelBegin(9), SG.outLabelBegin(9) + 4),
                      (std::vector<index>{5, 20, 40, 60}))
                << "directed=" << directed;

            // An arc that does not exist has no label, and nor does the reverse of a directed one.
            EXPECT_EQ(SG.edgeLabel(9, 7), none) << "directed=" << directed;
            EXPECT_EQ(SG.edgeLabel(6, 9), directed ? none : index{60}) << "directed=" << directed;

            // Nothing was thrown away, so no algorithm has anything to refuse.
            EXPECT_FALSE(SG.collapsedLabelledEdges()) << "directed=" << directed;
        }
    }

    // Then the compaction. It rewrites the slices in place and would corrupt the pairing just as
    // thoroughly as the sort if it moved a head without its label. Parallel arcs here all agree on
    // their label, so collapsing them is lossless and the survivor must keep it; self-loops are
    // dropped outright, whatever they carry.
    for (bool directed : {false, true}) {
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

        const IsomorphismDetails::SearchGraph SG(labelled.G, /* buildMatrix = */ true,
                                                 labelled.edgeLabels);

        expectSlicesAreSimpleNeighborhoods(labelled.G, SG);
        expectLabelsPairWithHeads(labelled.G, labelled.edgeLabels, SG);

        EXPECT_EQ(std::vector<node>(SG.outBegin(4), SG.outEnd(4)), (std::vector<node>{1, 3}))
            << "directed=" << directed;
        EXPECT_EQ(SG.edgeLabel(4, 1), 10u) << "directed=" << directed;
        EXPECT_EQ(SG.edgeLabel(4, 3), 30u) << "directed=" << directed;

        // Self-loops leave nothing behind, so their labels cannot disagree about anything.
        EXPECT_EQ(SG.outDegree(2), 0u) << "directed=" << directed;
        EXPECT_EQ(SG.edgeLabel(2, 2), none) << "directed=" << directed;
        EXPECT_FALSE(SG.collapsedLabelledEdges()) << "directed=" << directed;
    }

    // Last the two directions of a mutual pair. Those are two edges with two ids sitting in two
    // different slices, so they carry unrelated labels. Storing a label per node pair rather than
    // per arc collides them, and the assertions below say so immediately.
    const IsomorphismTest::LabelledGraph mutual = IsomorphismTest::labelledGraphOf(
        3, {{0, 1, 7}, {1, 0, 8}, {1, 2, 9}}, /* directed = */ true);

    const IsomorphismDetails::SearchGraph SM(mutual.G, /* buildMatrix = */ true, mutual.edgeLabels);

    EXPECT_EQ(SM.edgeLabel(0, 1), 7u);
    EXPECT_EQ(SM.edgeLabel(1, 0), 8u);
    EXPECT_EQ(SM.edgeLabel(1, 2), 9u);
    EXPECT_EQ(SM.edgeLabel(2, 1), none);
    EXPECT_FALSE(SM.collapsedLabelledEdges());

    expectLabelsPairWithHeads(mutual.G, mutual.edgeLabels, SM);

    // The in-slices are a second copy of the same arcs, so they have to agree arc for arc rather
    // than merely hold the right multiset. Node 1 is reached by 0 -> 1, which carries 7.
    ASSERT_EQ(SM.inDegree(1), 1u);
    EXPECT_EQ(*SM.inBegin(1), 0u);
    EXPECT_EQ(*SM.inLabelBegin(1), 7u);

    // ... and node 0 by 1 -> 0, which carries 8. Keying by node pair would give both the same.
    ASSERT_EQ(SM.inDegree(0), 1u);
    EXPECT_EQ(*SM.inBegin(0), 1u);
    EXPECT_EQ(*SM.inLabelBegin(0), 8u);
}

TEST_F(SearchGraphGTest, testCollapsedLabelledEdges) {

    // The flag every algorithm checks before it starts searching. It has to be exact in both
    // directions: never firing means a differently-labelled parallel edge is silently dropped and
    // the search answers a question nobody asked, while firing too eagerly refuses inputs that
    // collapse perfectly well.
    for (bool directed : {false, true}) {
        const IsomorphismTest::LabelledGraph plain =
            IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}}, directed);
        const IsomorphismTest::LabelledGraph agreeing = IsomorphismTest::labelledGraphOf(
            4, {{0, 1, 1}, {0, 1, 1}, {1, 2, 2}, {2, 3, 1}}, directed);
        const IsomorphismTest::LabelledGraph disagreeing = IsomorphismTest::labelledGraphOf(
            4, {{0, 1, 1}, {0, 1, 4}, {1, 2, 2}, {2, 3, 1}}, directed);

        EXPECT_FALSE(IsomorphismDetails::SearchGraph(plain.G, false, plain.edgeLabels)
                         .collapsedLabelledEdges())
            << "directed=" << directed << " - nothing was collapsed at all";

        EXPECT_FALSE(IsomorphismDetails::SearchGraph(agreeing.G, false, agreeing.edgeLabels)
                         .collapsedLabelledEdges())
            << "directed=" << directed << " - collapsing equal labels is lossless";

        EXPECT_TRUE(IsomorphismDetails::SearchGraph(disagreeing.G, false, disagreeing.edgeLabels)
                        .collapsedLabelledEdges())
            << "directed=" << directed << " - a label was thrown away";

        // An unlabelled snapshot of the very same graph has no labels to lose, so it must not
        // refuse anything. This is what keeps the flag from leaking into unlabelled searches.
        const IsomorphismDetails::SearchGraph unlabelled(disagreeing.G, false);
        EXPECT_FALSE(unlabelled.hasEdgeLabels()) << "directed=" << directed;
        EXPECT_FALSE(unlabelled.collapsedLabelledEdges()) << "directed=" << directed;
        EXPECT_EQ(unlabelled.edgeLabel(0, 1), none) << "directed=" << directed;
        EXPECT_EQ(unlabelled.outLabelBegin(0), nullptr) << "directed=" << directed;
    }
}

TEST_F(SearchGraphGTest, testEdgeLabelsNeedEdgeIds) {

    // Labels are indexed by edge id, so a graph without ids has no index space for them and the
    // scatter would read edgeLabels[none]. Saying so here beats crashing somewhere downstream.
    Graph G = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}});
    ASSERT_FALSE(G.hasEdgeIds());
    EXPECT_THROW(IsomorphismDetails::SearchGraph(G, false, std::vector<index>{1, 2}),
                 std::runtime_error);

    G.indexEdges();
    EXPECT_THROW(IsomorphismDetails::SearchGraph(G, false, std::vector<index>{1}),
                 std::runtime_error);
    EXPECT_NO_THROW(IsomorphismDetails::SearchGraph(G, false, std::vector<index>{1, 2}));
}

} // namespace NetworKit
