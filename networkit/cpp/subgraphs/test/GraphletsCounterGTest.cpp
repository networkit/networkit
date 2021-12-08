#include <array>

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/subgraphs/GraphletsCounter.hpp>

namespace NetworKit {

class GraphletsCounterGTest: public testing::Test {
};

static constexpr std::array<count, 5> nbGraphlets{{0, 0, 2, 4, 11}};

static Graph createUndirectedCompleteGraph(count V) {
    Graph K(
        V,
        false,  // weighted
        false   // directed
    );
    K.forNodePairs([&K](node u, node v) {
        K.addEdge(u, v);
    });
    return K;
}

static Graph createUndirectedCycle(count V) {
    Graph C(
        V,
        false,
        false
    );
    C.forNodes(
        [&C, V](node u) {
            C.addEdge(u, (u+1) % V);
        }
    );
    return C;
}

inline static count sum(const std::vector<count>& v) {
    count ret{0};
    for(count value : v)
        ret += value;
    return ret;
}

inline static count binom2(count n) {
    return (n*(n-1)) / 2;
}

inline static count binom3(count n) {
    return (n*(n-1)*(n-2)) / 6;
}

inline static count binom4(count n) {
    return (n*(n-1)*(n-2)*(n-3)) / 24;
}

inline static std::vector<count> expected3CountsForCompleteGraph(count V) {
    /* Suppose V > 2 */
    // Every subgraph of size 3 of K_n is a triangle
    // and there are n choose 3 such subgraphs
    return {binom3(V), 0, 0, 0};
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
    if(V == 3)
        return {1, 0, 0, 0};
    else if(V == 4)
        return {0, 4, 0, 0};
    else if(V == 5)
        return {0, 5, 5, 0};
    else  // V >= 6
        return {0, V, V*(V-4), binom3(V) - V*(V-3)};
}

inline static void matchExpected(const std::vector<count>& computed,
                                 const std::vector<count>& expected,
                                 count expected_size, count expected_sum) {
    ASSERT_EQ(expected.size(), expected_size) << "expected.size() == " << expected.size();
    ASSERT_EQ(sum(expected), expected_sum) << "sum(expected) == " << sum(expected);
    EXPECT_EQ(computed.size(), expected_size);
    for(index i{0}; i < static_cast<index>(expected_size); ++i)
        EXPECT_EQ(computed.at(i), expected.at(i));
}

TEST_F(GraphletsCounterGTest, testThrowingErrorForDirectedGraph) {
    Graph G(
        5,
        false,  // weighted
        true    // directed
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
    for(count V{3}; V <= 10; ++V) {
        Graph K{createUndirectedCompleteGraph(V)};
        GraphletsCounter counter(K, 3);
        counter.run();
        auto size3Graphlets{counter.getGraphletsCounts()};
        auto expected{expected3CountsForCompleteGraph(V)};
        matchExpected(
            size3Graphlets,
            expected3CountsForCompleteGraph(V),
            nbGraphlets.at(3),
            binom3(V)
        );
    }
}

TEST_F(GraphletsCounterGTest, test3GraphletsOnCycles) {
    for(count V{3}; V <= 10; ++V) {
        Graph C{createUndirectedCycle(V)};
        GraphletsCounter counter(C, 3);
        counter.run();
        auto size3Graphlets{counter.getGraphletsCounts()};
        auto expected{expected3CountsForCycle(V)};
        matchExpected(
            size3Graphlets,
            expected3CountsForCycle(V),
            nbGraphlets.at(3),
            binom3(V)
        );
    }
}

}  // namespacec NetworKit
