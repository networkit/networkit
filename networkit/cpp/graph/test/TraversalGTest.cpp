#include <gtest/gtest.h>

#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class TraversalGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    Graph generateRandomWeights(const Graph &G) const;
    bool weighted() const noexcept;
    bool directed() const noexcept;
};

INSTANTIATE_TEST_CASE_P(
    InstantiationName, TraversalGTest,
    testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                    std::make_pair(false, true),
                    std::make_pair(true, true)), ); // comma required for variadic macro

Graph TraversalGTest::generateRandomWeights(const Graph &G) const {
    Graph Gw(G, true, G.isDirected());
    Gw.forEdges([&](node u, node v) { Gw.setWeight(u, v, Aux::Random::probability()); });
    return Gw;
}

bool TraversalGTest::weighted() const noexcept { return GetParam().first; }

bool TraversalGTest::directed() const noexcept { return GetParam().second; }

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
            Traversal::BFSEdgesFrom(
                G, randNodes[i - 1], [&](node u, node v, edgeweight, edgeid) {
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

} // namespace NetworKit
