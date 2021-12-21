#include <array>
#include <numeric>

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/subgraphs/GraphletsCounter.hpp>

namespace NetworKit {

class GraphletsCounterGTest : public testing::Test {};

static constexpr std::array<count, 5> nbGraphlets{{0, 0, 2, 4, 11}};

static Graph createUndirectedCompleteGraph(count V) {
    Graph K(V,
            false, // weighted
            false  // directed
    );
    K.forNodePairs([&K](node u, node v) { K.addEdge(u, v); });
    return K;
}

static Graph createUndirectedCycle(count V) {
    Graph C(V, false, false);
    C.forNodes([&C, V](node u) { C.addEdge(u, (u + 1) % V); });
    return C;
}

/*
 * Create a graph G such that G = C_V \boxtimes C_V
 */
static Graph createUndirectedStrongProductOfCycles(count V) {
    Graph G(V * V, false, false);
    for (node u{0}; u < V; ++u) {
        node v{(u + 1) % V};
        for (node w{0}; w < V; ++w) {
            node x{(w + 1) % V};
            node uw{u * V + w}, vw{v * V + w}, ux{u * V + x}, vx{v * V + x};
            G.addEdge(uw, ux); // vertical
            G.addEdge(uw, vw); // horizontal
            G.addEdge(ux, vw); // main diagonal
            G.addEdge(uw, vx); // secondary diagonal
        }
    }
    return G;
}

inline static count sum(const std::vector<count> &v) {
    return std::accumulate(v.begin(), v.end(), static_cast<count>(0));
}

inline static count binom2(count n) {
    return (n * (n - 1)) / 2;
}

inline static count binom3(count n) {
    return (n * (n - 1) * (n - 2)) / 6;
}

inline static count binom4(count n) {
    return (n * (n - 1) * (n - 2) * (n - 3)) / 24;
}

inline static std::vector<count> expected3CountsForCompleteGraph(count V) {
    /* Suppose V > 2 */
    // Every subgraph of size 3 of K_n is a triangle
    // and there are n choose 3 such subgraphs
    return {binom3(V), 0, 0, 0};
}

inline static std::vector<count> expected4CountsForCompleteGraph(count V) {
    /* Suppose V > 3 */
    // Every subgraph of size 3 of K_n is a triangle
    // and there are n choose 3 such subgraphs
    // Every subgraph of size 4 of K_n is a 4-clique
    // and there are n choose 4 such graphs
    return {
        binom3(V), 0, 0, 0,                     // 3-graphlets
        binom4(V), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 // 4-graphlets
    };
}

inline static std::vector<count> expected3CountsForCycle(count V) {
    /* Suppose V > 2 */
    // A size 3 subgraph of a cycle (of > 5 vertices) is either:
    // (i)   a star;
    // (ii)  a disconnected graph with a single edge.
    // (iii) an empty graph
    // There are:
    // (i)   exactly V 2-stars (bijection between V and 2-stars
    //       by mapping every vertex v to the subgraph induced by
    //       {v, v+1, v+2} (mod V))
    // (ii)  exactly V(V-4) disconnected subgraphs with a single
    //       edge (choose any vertex v, add the vertex v+1 mod V and
    //       pick any vertex not in the range [v-1:v+2] mod V)
    // (iii) exactly V choose 3 - V(V-3) empty graphs (the only remaining ones)
    if (V == 3)
        return {1, 0, 0, 0};
    else
        return {0, V, V * (V - 4), binom3(V) - V * (V - 3)};
}

inline static std::vector<count> expected4CountsForCycle(count V) {
    /* Suppose V > 3 */
    // A size 4 subgraph of a cycle can never be either a
    // (i)   4-clique;
    // (ii)  4-chordal cycle
    // (iii) 4-tailed triangle
    // (iv)  3-star
    // because a cycle has ne vertex with a degree > 2.
    // Neither can it be a triangle with a disconnected vertex
    // (the complementary of a 3-star) because a subgraph of a path is
    // a path iff it is the whole graph
    // It can be:
    // (i)   a 4-cycle iff V == 4
    // (ii)  a 4-path if V > 4
    // (iii) a 2-star with an a disconnected vertex if V > 5
    // (iv)  a 2-disconnected edges graph if V > 5
    // (v)   a single edge with two disconnected vertices if V > 6
    // (vi)  the empty 4-graph if V > 7
    //
    // If V > 4, there are exactly V 4-path subgraphs in C_V:
    // pick any vertex v and consider the path {v, v+1, v+2, v+3} mod V
    //
    // If V > 5, there are exactly V(V-5) 2-stars with an added vertex:
    // pick any vertex v and consider the path (2-star) {v, v+1, v+2}
    // and then pick any of the V-5 remaining vertices by ignoring
    // {v-1, ..., v+3}. This last vertex is not connected any of the first 3.
    // There are also exactly V(V-5)/2 2-disconnected edges graphs:
    // pick any vertex v and consider the edge {v, v+1}. Then pick any
    // unselected vertex w, except {v-1, ..., v+2}. The subgraph induced by
    // {v, v+1, w, w+1} contains only 2 edges: {v, v+1} and {w, w+1}.
    // Indeed {v+1, w} is not an edge because w != v+2 and w != v
    // and {w+1, v} is not an edge because w != v-2 and w != v.
    // This accounts to V choices for V and V-5 choices for w
    // that still needs to be divided by 2 because first choosing v
    // then w or first w then v gives the same subgraph.
    //
    // If > 6, there are exactly V(V-5)(V-6)/2 subgraphs such that E = 1:
    // pick any vertex v and consider the edge {v, v+1}. The subgraph
    // induced by the remaining vertices when excluding {v-1, ..., v+2}
    // is a path of length V-5 in which we still need to pick to unconnected vertices.
    // Two vertices are unconnected in a graph iff they are in the omplementary
    // of the graph, i.e. there are V(V-1)/2 - E pairs of unconnected vertices in any
    // graph and in particular there are V(V-1)/2 - (V-1) = (V-1)(V-2)/2
    // unconnected pairs of vertices in a path of length V.
    // Therefore the total number of subgraphs such that E=1 is V multiplied
    // by ((V-4)-1) choose 2.
    //
    // Finally, a subgraph is the empty graph if it is none of the above.
    std::vector<count> ret(15, 0);
    auto size3Graphlets{expected3CountsForCycle(V)};
    std::copy(size3Graphlets.begin(), size3Graphlets.end(), ret.begin());
    if (V == 4)
        ret[7] = 1; // 4-cycle
    if (V > 4)
        ret[9] = V; // 4-path
    if (V > 5) {
        ret[11] = V * (V - 5);       // 2-star + 1 vertex
        ret[12] = (V * (V - 5)) / 2; // 2 disconnected edges
    }
    if (V > 6)
        ret[13] = V * binom2(V - 5);
    if (V > 7)
        ret[14] = binom3(V) + binom4(V) - sum(ret);
    return ret;
}

inline static std::vector<count> expected3CountsForStrongProductOfCycles(count V) {
    return {4 * V * V, 16 * V * V, 4 * V * V * (V * V - 13),
            binom3(V * V) - 4 * V * V * (V * V - 8)};
}

inline static std::vector<count> expected4CountsForStrongProductOfCycles(count V) {
    /* Suppose V > 4 */
    std::vector<count> ret(15, 0);
    auto size3Graphlets{expected3CountsForStrongProductOfCycles(V)};
    std::copy(size3Graphlets.begin(), size3Graphlets.end(), ret.begin());
    ret[4] = V * V;
    ret[5] = 8 * V * V;
    ret[6] = 28 * V * V;
    ret[7] = V * V;
    ret[8] = 8 * V * V;
    ret[9] = 64 * V * V;
    ret[10] = 4 * V * V * (V * V - 15);
    ret[11] = 4 * V * V * (4 * V * V - 69);
    ret[12] = V * V * (8 * V * V - 143);
    ret[13] = 2 * V * V * (V * V * V * V - 35 * V * V + 326);
    ret[14] = binom3(V * V) + binom4(V * V) - sum(ret);
    return ret;
}

inline static void matchExpected(const std::vector<count> &computed,
                                 const std::vector<count> &expected, count expected_size,
                                 count expected_sum) {
    ASSERT_EQ(expected.size(), expected_size);
    ASSERT_EQ(sum(expected), expected_sum);
    EXPECT_EQ(computed.size(), expected_size);
    for (index i{0}; i < static_cast<index>(expected_size); ++i)
        EXPECT_EQ(computed.at(i), expected.at(i));
}

TEST_F(GraphletsCounterGTest, testThrowingErrorForDirectedGraph) {
    Graph G(5,
            false, // weighted
            true   // directed
    );
    GraphletsCounter counter(G, 2);
    EXPECT_THROW(counter.run(), std::runtime_error);
}

TEST_F(GraphletsCounterGTest, testThrowingErrorForWrongK) {
    Graph G{createUndirectedCompleteGraph(10)};
    GraphletsCounter counter(G, 0u);
    EXPECT_THROW(counter.run(), std::runtime_error);
}

TEST_F(GraphletsCounterGTest, test3GraphletsOnCompleteGraphs) {
    for (count V{3}; V <= 20; ++V) {
        Graph K{createUndirectedCompleteGraph(V)};
        GraphletsCounter counter(K, 3);
        counter.run();
        matchExpected(counter.getGraphletsCounts(), expected3CountsForCompleteGraph(V),
                      nbGraphlets.at(3), binom3(V));
    }
}

TEST_F(GraphletsCounterGTest, test3GraphletsOnCycles) {
    for (count V{3}; V <= 20; ++V) {
        Graph C{createUndirectedCycle(V)};
        GraphletsCounter counter(C, 3);
        counter.run();
        matchExpected(counter.getGraphletsCounts(), expected3CountsForCycle(V), nbGraphlets.at(3),
                      binom3(V));
    }
}

TEST_F(GraphletsCounterGTest, test4GraphletsOnCompleteGraphs) {
    for (count V{4}; V <= 20; ++V) {
        Graph K{createUndirectedCompleteGraph(V)};
        GraphletsCounter counter(K, 4);
        counter.run();
        matchExpected(counter.getGraphletsCounts(), expected4CountsForCompleteGraph(V),
                      nbGraphlets.at(3) + nbGraphlets.at(4), binom3(V) + binom4(V));
    }
}

TEST_F(GraphletsCounterGTest, test4GraphletsOnCycles) {
    for (count V{4}; V <= 20; ++V) {
        Graph C{createUndirectedCycle(V)};
        GraphletsCounter counter(C, 4);
        counter.run();
        matchExpected(counter.getGraphletsCounts(), expected4CountsForCycle(V),
                      nbGraphlets.at(3) + nbGraphlets.at(4), binom3(V) + binom4(V));
    }
}

TEST_F(GraphletsCounterGTest, test3GraphletsOnStrongProductsOfCycles) {
    for (count V{4}; V <= 20; ++V) {
        Graph Csquared{createUndirectedStrongProductOfCycles(V)};
        GraphletsCounter counter(Csquared, 3);
        counter.run();
        matchExpected(counter.getGraphletsCounts(), expected3CountsForStrongProductOfCycles(V),
                      nbGraphlets.at(3), binom3(V * V));
    }
}

TEST_F(GraphletsCounterGTest, test4GraphletsOnStrongProductsOfCycles) {
    for (count V{5}; V <= 20; ++V) {
        Graph Csquared{createUndirectedStrongProductOfCycles(V)};
        GraphletsCounter counter(Csquared, 4);
        counter.run();
        matchExpected(counter.getGraphletsCounts(), expected4CountsForStrongProductOfCycles(V),
                      nbGraphlets.at(3) + nbGraphlets.at(4), binom3(V * V) + binom4(V * V));
    }
}

} // namespace NetworKit
