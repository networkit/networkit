#include <cstdint>
#include <queue>
#include <stack>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

template <class GraphT_>
struct BFSTraversalConfig {
    using GraphT = GraphT_;
};

template <class TestT>
class GenericBFSTraversalGTest : public testing::Test {
public:
    using GraphT = typename TestT::GraphT;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    GraphT unweightedGraph() const {
        GraphT G(7, false);
        G.addEdge(NodeT{0}, NodeT{1});
        G.addEdge(NodeT{0}, NodeT{2});
        G.addEdge(NodeT{1}, NodeT{3});
        G.addEdge(NodeT{2}, NodeT{4});
        G.addEdge(NodeT{3}, NodeT{5});
        G.addEdge(NodeT{4}, NodeT{6});
        return G;
    }

    GraphT weightedGraph() const {
        GraphT G(7, true);
        G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{10});
        G.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{20});
        G.addEdge(NodeT{1}, NodeT{3}, EdgeWeightT{30});
        G.addEdge(NodeT{2}, NodeT{4}, EdgeWeightT{40});
        G.addEdge(NodeT{3}, NodeT{5}, EdgeWeightT{50});
        G.addEdge(NodeT{4}, NodeT{6}, EdgeWeightT{60});
        G.addEdge(NodeT{2}, NodeT{3}, EdgeWeightT{70});
        return G;
    }
};

using BFSTraversalTestTypes =
    ::testing::Types<BFSTraversalConfig<Graph>, BFSTraversalConfig<AdjListGraph<uint32_t, float>>,
                     BFSTraversalConfig<AdjListGraph<int, int>>>;

TYPED_TEST_SUITE(GenericBFSTraversalGTest, BFSTraversalTestTypes);

TYPED_TEST(GenericBFSTraversalGTest, testBFSfromTypedGraphs) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->unweightedGraph();

    std::vector<NodeT> sourceSequence;
    std::vector<count> sourceDistances;
    Traversal::BFSfrom(G, NodeT{0}, [&](NodeT u, count dist) {
        sourceSequence.push_back(u);
        sourceDistances.push_back(dist);
    });

    EXPECT_THAT(sourceSequence, testing::ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{3},
                                                     NodeT{4}, NodeT{5}, NodeT{6}));
    EXPECT_THAT(sourceDistances, testing::ElementsAre(0, 1, 1, 2, 2, 3, 3));

    const std::vector<NodeT> sources{NodeT{0}, NodeT{6}};
    std::vector<NodeT> rangeSequence;
    Traversal::BFSfrom(G, sources.begin(), sources.end(),
                       [&](NodeT u) { rangeSequence.push_back(u); });

    EXPECT_THAT(rangeSequence, testing::ElementsAre(NodeT{0}, NodeT{6}, NodeT{1}, NodeT{2},
                                                    NodeT{4}, NodeT{3}, NodeT{5}));
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSEdgesFromTypedGraphs) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    const auto G = this->weightedGraph();
    std::vector<std::tuple<NodeT, NodeT, EdgeWeightT>> edgeSequence;

    Traversal::BFSEdgesFrom(G, NodeT{0}, [&](NodeT u, NodeT v, EdgeWeightT w, edgeid) {
        edgeSequence.emplace_back(u, v, w);
    });

    EXPECT_THAT(edgeSequence,
                testing::ElementsAre(std::make_tuple(NodeT{0}, NodeT{1}, EdgeWeightT{10}),
                                     std::make_tuple(NodeT{0}, NodeT{2}, EdgeWeightT{20}),
                                     std::make_tuple(NodeT{1}, NodeT{3}, EdgeWeightT{30}),
                                     std::make_tuple(NodeT{2}, NodeT{4}, EdgeWeightT{40}),
                                     std::make_tuple(NodeT{3}, NodeT{5}, EdgeWeightT{50}),
                                     std::make_tuple(NodeT{4}, NodeT{6}, EdgeWeightT{60})));
}

class TraversalGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool weighted() const noexcept;
    bool directed() const noexcept;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, TraversalGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

bool TraversalGTest::weighted() const noexcept {
    return GetParam().first;
}

bool TraversalGTest::directed() const noexcept {
    return GetParam().second;
}

TEST_P(TraversalGTest, testBFSfrom) {
    constexpr count n = 200;
    constexpr double p = 0.15;
    std::vector<unsigned char> visited(n);
    std::vector<node> sequence;
    sequence.reserve(n);
    std::vector<std::pair<node, node>> edgeSequence;

    auto doBFS = [&](const Graph &G, const std::vector<node> &sources) {
        std::fill(visited.begin(), visited.end(), 0);
        sequence.clear();
        edgeSequence.clear();
        std::queue<node> q;

        for (node source : sources) {
            q.push(source);
            visited[source] = 1;
        }

        do {
            node u = q.front();
            q.pop();
            sequence.push_back(u);
            G.forNeighborsOf(u, [&](node v) {
                if (!visited[v]) {
                    q.push(v);
                    visited[v] = 1;
                    edgeSequence.push_back({u, v});
                }
            });
        } while (!q.empty());
    };

    std::vector<node> randNodes;
    for (node u = 0; u < n; ++u) {
        randNodes.push_back(u);
    }

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        std::shuffle(randNodes.begin(), randNodes.end(), Aux::Random::getURNG());
        const auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        for (count i = 1; i <= n; ++i) {
            std::vector<node> sources(randNodes.begin(), randNodes.begin() + i);
            doBFS(G, sources);
            count curNode = 0;
            Traversal::BFSfrom(G, sources.begin(), sources.end(),
                               [&](node u) { EXPECT_EQ(sequence[curNode++], u); });

            sources.clear();
            sources.push_back(randNodes[i - 1]);
            doBFS(G, sources);
            curNode = 0;
            Traversal::BFSEdgesFrom(G, randNodes[i - 1], [&](node u, node v, edgeweight, edgeid) {
                EXPECT_EQ(edgeSequence[curNode++], std::make_pair(u, v));
            });
        }
    }
}

TEST_P(TraversalGTest, testDFSfrom) {
    constexpr count n = 200;
    constexpr double p = 0.15;
    std::vector<unsigned char> visited(n);
    std::vector<node> sequence;
    sequence.reserve(n);
    std::vector<std::pair<node, node>> edgeSequence;

    auto doDFS = [&](const Graph &G, node source) {
        sequence.clear();
        edgeSequence.clear();
        std::fill(visited.begin(), visited.end(), 0);
        visited[source] = 1;
        std::stack<node> s;
        s.push(source);

        do {
            node u = s.top();
            s.pop();
            sequence.push_back(u);
            G.forNeighborsOf(u, [&](node v) {
                if (!visited[v]) {
                    s.push(v);
                    visited[v] = 1;
                    edgeSequence.push_back({u, v});
                }
            });
        } while (!s.empty());
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        const auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        G.forNodes([&](node source) {
            doDFS(G, source);
            count curNode = 0;
            Traversal::DFSfrom(G, source, [&](node u) { EXPECT_EQ(sequence[curNode++], u); });
            curNode = 0;
            Traversal::DFSEdgesFrom(G, source, [&](node u, node v, edgeweight, edgeid) {
                EXPECT_EQ(edgeSequence[curNode++], std::make_pair(u, v));
            });
        });
    }
}

TEST_P(TraversalGTest, testDijkstraFrom) {
    constexpr count n = 200;
    constexpr double p = 0.15;

    std::vector<node> randNodes, nodes;
    for (node u = 0; u < n; ++u) {
        randNodes.push_back(u);
        nodes.push_back(u);
    }

    auto dijkstra = [&](const Graph &G, const std::vector<node> &sources) {
        std::vector<edgeweight> distance(G.upperNodeIdBound(),
                                         std::numeric_limits<edgeweight>::max());
        tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<edgeweight>> heap{distance};
        for (const auto u : sources) {
            distance[u] = 0;
            heap.push(u);
        }

        do {
            const auto u = heap.extract_top();
            G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
                if (distance[v] > distance[u] + w) {
                    distance[v] = distance[u] + w;
                    heap.update(v);
                }
            });
        } while (!heap.empty());

        return distance;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        if (weighted()) {
            GraphTools::randomizeWeights(G);
        }

        G.forNodes([&](const node u) {
            const auto distance = dijkstra(G, {u});
            std::sort(nodes.begin(), nodes.end(), [&distance](const node x, const node y) {
                return distance[x] < distance[y];
            });
            for (node i = 0; i < n - 1; ++i)
                assert(distance[nodes[i]] <= distance[nodes[i + 1]]);

            index i = 0;
            Traversal::DijkstraFrom(G, u, [&](const node u, const edgeweight w) {
                if (u != nodes[i]) {
                    EXPECT_DOUBLE_EQ(distance[u], distance[nodes[i]]);
                }
                ++i;
                EXPECT_DOUBLE_EQ(w, distance[u]);
            });
        });

        std::shuffle(randNodes.begin(), randNodes.end(), Aux::Random::getURNG());
        for (count i = 1; i <= n; ++i) {
            const auto distance =
                dijkstra(G, std::vector<node>(randNodes.begin(), randNodes.begin() + i));
            Traversal::DijkstraFrom(G, randNodes.begin(), randNodes.begin() + i,
                                    [&distance](const node u, const edgeweight d) {
                                        EXPECT_DOUBLE_EQ(d, distance[u]);
                                    });
        }
    }
}

} // namespace NetworKit
