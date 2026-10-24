#include <algorithm>
#include <cstdint>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISGraphWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

// Prints a summary in failure messages instead of the raw bytes of the graph. Must be in the
// namespace of the graph to be found by gtest, static to stay local to this file.
template <class NodeT, class EdgeWeightT>
static void PrintTo(const AdjListGraph<NodeT, EdgeWeightT> &G, std::ostream *os) {
    *os << (G.isDirected() ? "directed" : "undirected") << ' '
        << (G.isWeighted() ? "weighted" : "unweighted") << " graph with " << G.numberOfNodes()
        << " nodes and " << G.numberOfEdges() << " edges";
}

namespace {

using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

// Matches a graph of the same type with the same nodes, edges, and edge weights as `expected`.
// Explains the first difference.
MATCHER_P(IsSameGraphAs, expected, "has the same nodes, edges, and edge weights") {
    using GraphT = std::remove_cvref_t<decltype(arg)>;
    static_assert(std::is_same_v<GraphT, std::remove_cvref_t<decltype(expected)>>,
                  "IsSameGraphAs compares graphs of the same type only");
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    const auto explainCount = [result_listener](const char *what, count actual,
                                                count wanted) -> bool {
        *result_listener << what << " is " << actual << " instead of " << wanted;
        return false;
    };
    if (arg.numberOfNodes() != expected.numberOfNodes()) {
        return explainCount("number of nodes", arg.numberOfNodes(), expected.numberOfNodes());
    }
    if (arg.upperNodeIdBound() != expected.upperNodeIdBound()) {
        return explainCount("upper node id bound", static_cast<count>(arg.upperNodeIdBound()),
                            static_cast<count>(expected.upperNodeIdBound()));
    }
    if (arg.numberOfEdges() != expected.numberOfEdges()) {
        return explainCount("number of edges", arg.numberOfEdges(), expected.numberOfEdges());
    }
    if (arg.isWeighted() != expected.isWeighted()) {
        *result_listener << "is " << (arg.isWeighted() ? "weighted" : "unweighted");
        return false;
    }
    if (arg.isDirected() != expected.isDirected()) {
        *result_listener << "is " << (arg.isDirected() ? "directed" : "undirected");
        return false;
    }

    // Same number of edges, so it suffices to check that every expected edge is in arg.
    const auto isInArg = [&](const WeightedEdgeT<NodeT, EdgeWeightT> &edge) -> bool {
        if (!arg.hasEdge(edge.u, edge.v)) {
            *result_listener << "has no edge " << edge.u << "-" << edge.v;
            return false;
        }
        if (arg.weight(edge.u, edge.v) != edge.weight) {
            // The cast prints 8-bit integral weights as numbers instead of characters.
            using PrintedWeightT = std::common_type_t<EdgeWeightT, int>;
            *result_listener << "edge " << edge.u << "-" << edge.v << " has weight "
                             << static_cast<PrintedWeightT>(arg.weight(edge.u, edge.v))
                             << " instead of " << static_cast<PrintedWeightT>(edge.weight);
            return false;
        }
        return true;
    };
    return std::ranges::all_of(expected.edgeWeightRange(), isInArg);
}

template <class GraphT_>
class METISGraphIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

protected:
    // Returns a random graph with weights representable by all tested edge weight types, or without
    // weights. It is sparse enough to have isolated nodes, whose METIS adjacency lines are empty.
    static GraphT randomGraph(bool weighted) {
        constexpr count n = 50;
        Aux::Random::setSeed(42, false);
        const Graph edges = ErdosRenyiGenerator(n, 0.05).generate();
        GraphT G(n, weighted);
        edges.forEdges([&](node u, node v) -> void {
            G.addEdge(static_cast<NodeT>(u), static_cast<NodeT>(v),
                      static_cast<EdgeWeightT>(Aux::Random::integer(1, 100)));
        });
        return G;
    }

    // Returns the contents of a METIS file storing G. Does not use the METIS writer, so that the
    // read tests do not depend on it.
    static std::string toMETIS(const GraphT &G) {
        std::ostringstream metis;
        metis << "% A comment\n"
              << G.numberOfNodes() << ' ' << G.numberOfEdges() << ' ' << int{G.isWeighted()}
              << '\n';
        G.forNodes([&](NodeT u) -> void {
            G.forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT w, edgeid) -> void {
                metis << static_cast<count>(v) + 1 << ' ';
                if (G.isWeighted()) {
                    // The cast prints 8-bit integral weights as numbers instead of characters.
                    metis << static_cast<std::common_type_t<EdgeWeightT, int>>(w) << ' ';
                }
            });
            metis << '\n';
        });
        return metis.str();
    }
};

using METISGraphIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint32_t, uint8_t>>;

TYPED_TEST_SUITE(METISGraphIOGTest, METISGraphIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(METISGraphIOGTest, testReadUnweighted) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    const GraphT expected = TestFixture::randomGraph(false);
    const std::string path =
        directory.writeFile("unweighted.graph", TestFixture::toMETIS(expected));
    EXPECT_THAT(GenericMETISGraphReader<GraphT>{}.read(path), IsSameGraphAs(expected));
}

TYPED_TEST(METISGraphIOGTest, testReadIntegerWeights) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    const GraphT expected = TestFixture::randomGraph(true);
    const std::string path = directory.writeFile("weighted.graph", TestFixture::toMETIS(expected));
    EXPECT_THAT(GenericMETISGraphReader<GraphT>{}.read(path), IsSameGraphAs(expected));
}

TYPED_TEST(METISGraphIOGTest, testReadFractionalWeights) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    // A weighted triangle whose weights are exactly representable by float and double.
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("fractional_weights.graph", "3 3 1\n"
                                                                             "2 0.5 3 1.25\n"
                                                                             "1 0.5 3 2.75\n"
                                                                             "1 1.25 2 2.75\n");
    GenericMETISGraphReader<GraphT> reader;
    if constexpr (std::is_integral_v<EdgeWeightT>) {
        EXPECT_THAT([&]() -> void { reader.read(path); },
                    ThrowsMessage<std::runtime_error>(HasSubstr("METIS value 0.5 is not valid")));
    } else {
        GraphT expected(3, true);
        expected.addEdge(NodeT{0}, NodeT{1}, static_cast<EdgeWeightT>(0.5));
        expected.addEdge(NodeT{0}, NodeT{2}, static_cast<EdgeWeightT>(1.25));
        expected.addEdge(NodeT{1}, NodeT{2}, static_cast<EdgeWeightT>(2.75));
        EXPECT_THAT(reader.read(path), IsSameGraphAs(expected));
    }
}

TYPED_TEST(METISGraphIOGTest, testReadLargerGraph) {
    using GraphT = typename TestFixture::GraphT;
    const GraphT G = GenericMETISGraphReader<GraphT>{}.read("input/jazz.graph");
    EXPECT_FALSE(G.isWeighted());
    EXPECT_EQ(G.numberOfNodes(), count{198});
    EXPECT_EQ(G.numberOfEdges(), count{2742});
}

TYPED_TEST(METISGraphIOGTest, testWriteAndReadBack) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    const std::string path = directory.path("roundtrip.graph");

    GraphT G(6, true);
    G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{2});
    G.addEdge(NodeT{1}, NodeT{3}, EdgeWeightT{5});
    G.addEdge(NodeT{3}, NodeT{5}, EdgeWeightT{7});
    G.addEdge(NodeT{0}, NodeT{4}, EdgeWeightT{100});
    // Leave a gap in the node ids, the writer has to compact them.
    G.removeNode(NodeT{2});

    GenericMETISGraphWriter<GraphT>{}.write(G, path);

    // Nodes after the removed node 2 are shifted down by one.
    GraphT expected(5, true);
    expected.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{2});
    expected.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{5});
    expected.addEdge(NodeT{2}, NodeT{4}, EdgeWeightT{7});
    expected.addEdge(NodeT{0}, NodeT{3}, EdgeWeightT{100});
    EXPECT_THAT(GenericMETISGraphReader<GraphT>{}.read(path), IsSameGraphAs(expected));
}

TEST(METISGraphIOGTest, testReadThrowsIfNodesDoNotFitNodeType) {
    using Reader = GenericMETISGraphReader<AdjListGraph<uint16_t, double>>;

    constexpr count maxNodes = std::numeric_limits<uint16_t>::max();
    const TemporaryDirectory directory;

    // The largest node id maxNodes - 1 is still smaller than NullNodeId<uint16_t>.
    const std::string fits =
        directory.writeFile("max_nodes.graph", std::to_string(maxNodes) + " 0\n");
    EXPECT_EQ(Reader{}.read(fits).numberOfNodes(), maxNodes);

    const std::string tooMany =
        directory.writeFile("too_many_nodes.graph", std::to_string(maxNodes + 1) + " 0\n");
    EXPECT_THAT([&]() -> void { Reader{}.read(tooMany); },
                ThrowsMessage<std::runtime_error>(
                    HasSubstr("has " + std::to_string(maxNodes + 1)
                              + " nodes, but the node type of the graph can hold at most "
                              + std::to_string(maxNodes) + " nodes")));
}

} // namespace
} // namespace NetworKit
