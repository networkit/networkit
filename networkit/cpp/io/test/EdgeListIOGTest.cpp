#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/EdgeListWriter.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {

namespace {

using ::testing::AllOf;
using ::testing::HasSubstr;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

template <class GraphT_>
class EdgeListIOGTest : public ::testing::Test {
public:
    using GraphT = GraphT_;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

protected:
    // Checks that G has exactly the edges {0, 1}, {0, 2}, and {1, 2}, with the given weights.
    static void expectTriangle(const GraphT &G, bool weighted, EdgeWeightT w01, EdgeWeightT w02,
                               EdgeWeightT w12) {
        EXPECT_EQ(G.numberOfNodes(), count{3});
        EXPECT_EQ(G.numberOfEdges(), count{3});
        EXPECT_EQ(G.isWeighted(), weighted);
        EXPECT_TRUE(G.checkConsistency());
        EXPECT_EQ(G.weight(NodeT{0}, NodeT{1}), w01);
        EXPECT_EQ(G.weight(NodeT{0}, NodeT{2}), w02);
        EXPECT_EQ(G.weight(NodeT{1}, NodeT{2}), w12);
    }
};

using EdgeListIOTestTypes =
    ::testing::Types<Graph, AdjListGraph<uint32_t, float>, AdjListGraph<uint16_t, double>,
                     AdjListGraph<uint64_t, uint64_t>, AdjListGraph<int, int>,
                     AdjListGraph<int64_t, float>, AdjListGraph<uint32_t, uint8_t>>;

TYPED_TEST_SUITE(EdgeListIOGTest, EdgeListIOTestTypes, /*Comma needed for variadic macro.*/);

TYPED_TEST(EdgeListIOGTest, testReadUnweighted) {
    using GraphT = typename TestFixture::GraphT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("unweighted.edgelist", "# A comment\n"
                                                                        "1 2\n"
                                                                        "\n"
                                                                        "1 3\n"
                                                                        "3 2\n");
    const GraphT G = GenericEdgeListReader<GraphT>(' ', 1).read(path);
    TestFixture::expectTriangle(G, false, EdgeWeightT{1}, EdgeWeightT{1}, EdgeWeightT{1});
}

TYPED_TEST(EdgeListIOGTest, testReadWeighted) {
    using GraphT = typename TestFixture::GraphT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("weighted.edgelist", "0\t1\t2\n"
                                                                      "0\t2\t4\n"
                                                                      "1\t2\t3\n");
    const GraphT G = GenericEdgeListReader<GraphT>('\t', 0).read(path);
    TestFixture::expectTriangle(G, true, EdgeWeightT{2}, EdgeWeightT{4}, EdgeWeightT{3});
}

// Only the last line has a weight: the sweep finds it, and the other edges get weight 1.
TYPED_TEST(EdgeListIOGTest, testReadPartiallyWeighted) {
    using GraphT = typename TestFixture::GraphT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("partial.edgelist", "0 1\n"
                                                                     "0 2\n"
                                                                     "1 2 7\r\n");
    const GraphT G = GenericEdgeListReader<GraphT>(' ', 0).read(path);
    TestFixture::expectTriangle(G, true, EdgeWeightT{1}, EdgeWeightT{1}, EdgeWeightT{7});
}

TYPED_TEST(EdgeListIOGTest, testReadFractionalWeights) {
    using GraphT = typename TestFixture::GraphT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    // Weights exactly representable by float and double.
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("fractional.edgelist", "0 1 0.5\n"
                                                                        "0 2 1.25\n"
                                                                        "1 2 2.75\n");
    GenericEdgeListReader<GraphT> reader(' ', 0);
    if constexpr (std::is_integral_v<EdgeWeightT>) {
        EXPECT_THROW(reader.read(path), std::runtime_error);
    } else {
        TestFixture::expectTriangle(reader.read(path), true, EdgeWeightT{0.5}, EdgeWeightT{1.25},
                                    EdgeWeightT{2.75});
    }
}

TYPED_TEST(EdgeListIOGTest, testReadDirected) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("directed.edgelist", "0 1\n"
                                                                      "1 0\n"
                                                                      "1 2\n");
    const GraphT G = GenericEdgeListReader<GraphT>(' ', 0, "#", true, true).read(path);
    EXPECT_TRUE(G.isDirected());
    EXPECT_EQ(G.numberOfEdges(), count{3});
    EXPECT_TRUE(G.hasEdge(NodeT{1}, NodeT{2}));
    EXPECT_FALSE(G.hasEdge(NodeT{2}, NodeT{1}));
}

TYPED_TEST(EdgeListIOGTest, testReadNonContinuous) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("alphabet.edgelist", "# comment\n"
                                                                      "a b 3\n"
                                                                      "a c 2\n"
                                                                      "c d\n");
    GenericEdgeListReader<GraphT> reader(' ', 0, "#", false);
    const GraphT G = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), count{4});
    EXPECT_TRUE(G.isWeighted());
    EXPECT_EQ(G.weight(NodeT{0}, NodeT{1}), EdgeWeightT{3});
    EXPECT_EQ(G.weight(NodeT{0}, NodeT{2}), EdgeWeightT{2});
    EXPECT_EQ(G.weight(NodeT{2}, NodeT{3}), EdgeWeightT{1});
    EXPECT_EQ(reader.getNodeMap().at("d"), NodeT{3});
}

TYPED_TEST(EdgeListIOGTest, testWriteAndReadBack) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;
    const TemporaryDirectory directory;

    GraphT G(5, true);
    G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{2});
    G.addEdge(NodeT{1}, NodeT{3}, EdgeWeightT{5});
    G.addEdge(NodeT{3}, NodeT{4}, EdgeWeightT{100});
    G.addEdge(NodeT{0}, NodeT{4}, EdgeWeightT{7});

    for (const bool bothDirections : {false, true}) {
        for (const NodeT firstNode : {NodeT{0}, NodeT{1}}) {
            const std::string path = directory.path("roundtrip.edgelist");
            GenericEdgeListWriter<GraphT>(' ', firstNode, bothDirections).write(G, path);
            const GraphT read = GenericEdgeListReader<GraphT>(' ', firstNode).read(path);
            EXPECT_TRUE(read.checkConsistency());
            EXPECT_EQ(read.numberOfNodes(), count{5});
            EXPECT_EQ(read.numberOfEdges(), count{4});
            EXPECT_TRUE(read.isWeighted());
            G.forEdges([&](NodeT u, NodeT v, EdgeWeightT w, edgeid) -> void {
                EXPECT_TRUE(read.hasEdge(u, v));
                EXPECT_EQ(read.weight(u, v), w);
            });
        }
    }
}

TYPED_TEST(EdgeListIOGTest, testWriteUnweighted) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    const TemporaryDirectory directory;
    // Directed, so that the order of the nodes in a line is deterministic.
    GraphT G(3, false, true);
    G.addEdge(NodeT{0}, NodeT{1});
    G.addEdge(NodeT{1}, NodeT{2});

    const std::string path = directory.path("unweighted.edgelist");
    GenericEdgeListWriter<GraphT>('\t', 1).write(G, path);

    std::ifstream file(path);
    const std::string text{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
    EXPECT_EQ(text, "1\t2\n2\t3\n");
}

TYPED_TEST(EdgeListIOGTest, testReadThrowsIfNodeIdTooSmall) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("too_small.edgelist", "0 1\n");
    EXPECT_THAT([&]() -> void { GenericEdgeListReader<GraphT>(' ', 1).read(path); },
                ThrowsMessage<std::runtime_error>(HasSubstr("smaller than the first node")));
}

TYPED_TEST(EdgeListIOGTest, testReadThrowsIfLineIsMalformed) {
    using GraphT = typename TestFixture::GraphT;
    const TemporaryDirectory directory;
    GenericEdgeListReader<GraphT> reader(' ', 0);
    EXPECT_THROW(reader.read(directory.writeFile("letters.edgelist", "0 x\n")), std::runtime_error);
    EXPECT_THROW(reader.read(directory.writeFile("single.edgelist", "0\n")), std::runtime_error);
}

TEST(EdgeListIOGTest, testFractionalNodeIdThrows) {
    const TemporaryDirectory directory;
    EXPECT_THROW(EdgeListReader(' ', 0).read(directory.writeFile("fraction.edgelist", "0 1.5\n")),
                 std::runtime_error);
}

TEST(EdgeListIOGTest, testReadThrowsIfNodeIdDoesNotFitNodeType) {
    using Reader = GenericEdgeListReader<AdjListGraph<uint16_t, double>>;
    constexpr uint32_t maxNode = std::numeric_limits<uint16_t>::max();
    const TemporaryDirectory directory;

    // The largest node id maxNode - 1 is still smaller than NullNodeId<uint16_t>.
    const std::string fits =
        directory.writeFile("max_id.edgelist", "0 " + std::to_string(maxNode - 1) + "\n");
    EXPECT_EQ(Reader(' ', 0).read(fits).numberOfNodes(), count{maxNode});

    const std::string nullId =
        directory.writeFile("null_id.edgelist", "0 " + std::to_string(maxNode) + "\n");
    EXPECT_THAT([&]() -> void { Reader(' ', 0).read(nullId); },
                ThrowsMessage<std::runtime_error>(HasSubstr("does not fit")));

    const std::string tooLarge =
        directory.writeFile("large_id.edgelist", "0 " + std::to_string(maxNode + 1) + "\n");
    EXPECT_THROW(Reader(' ', 0).read(tooLarge), std::runtime_error);
}

TEST(EdgeListIOGTest, testWriteFloatingPointWeightsExactly) {
    using GraphT = AdjListGraph<uint32_t, double>;
    const TemporaryDirectory directory;
    const std::string path = directory.path("precise.edgelist");

    GraphT G(3, true);
    G.addEdge(0u, 1u, 1234567.0);
    G.addEdge(1u, 2u, 53.71);
    GenericEdgeListWriter<GraphT>(' ', 0).write(G, path);

    const GraphT read = GenericEdgeListReader<GraphT>(' ', 0).read(path);
    EXPECT_EQ(read.weight(0u, 1u), 1234567.0);
    EXPECT_EQ(read.weight(1u, 2u), 53.71);
}

} // namespace
} // namespace NetworKit
