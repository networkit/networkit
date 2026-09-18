/*
 * NetworkitBinaryGTest.cpp
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <random>
#include <string>
#include <variant>
#include <vector>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <networkit/io/NetworkitBinaryReader.hpp>
#include <networkit/io/NetworkitBinaryWriter.hpp>
#include <networkit/io/SNAPGraphReader.hpp>

namespace NetworKit {

class NetworkitBinaryGTest : public testing::Test {};

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryTiny01) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/tiny_01.graph");
    NetworkitBinaryWriter writer;

    writer.write(G, "output/binary_tiny01");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_tiny01");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    G.forNodes([&](node u) { G.forEdgesOf(u, [&](node v) { ASSERT_TRUE(G2.hasEdge(u, v)); }); });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryTiny01InMemory) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/tiny_01.graph");
    NetworkitBinaryWriter writer;

    std::vector<uint8_t> data = writer.writeToBuffer(G);
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.readFromBuffer(data);
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    G.forNodes([&](node u) { G.forEdgesOf(u, [&](node v) { ASSERT_TRUE(G2.hasEdge(u, v)); }); });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryAdjListGraphTemplateInMemory) {
    using SmallGraph = AdjListGraph<uint32_t, float>;

    SmallGraph G(6, true, true);
    G.addEdge(0, 1, 1.25f);
    G.addEdge(2, 4, -3.5f);
    G.addEdge(4, 2, 7.0f);
    G.addEdge(5, 5, 2.0f);
    G.indexEdges();

    NetworkitBinaryWriter writer(4, NetworkitBinaryWeights::AUTO_DETECT);
    const auto data = writer.writeToBuffer(G);

    ASSERT_EQ(0, std::memcmp(data.data(), "nkbg005", 8));
    uint64_t features;
    std::memcpy(&features, data.data() + 2 * sizeof(uint64_t), sizeof(uint64_t));
    EXPECT_EQ(nkbg::widthCode(sizeof(SmallGraph::NodeT)),
              (features & nkbg::NODE_TYPE_WIDTH_MASK) >> nkbg::NODE_TYPE_WIDTH_SHIFT);

    NetworkitBinaryReader reader;
    const SmallGraph G2 = reader.readFromBuffer<SmallGraph>(data);
    EXPECT_TRUE(G2.isDirected());
    EXPECT_TRUE(G2.isWeighted());
    EXPECT_TRUE(G2.hasEdgeIds());
    ASSERT_EQ(G.numberOfNodes(), G2.numberOfNodes());
    ASSERT_EQ(G.numberOfEdges(), G2.numberOfEdges());
    G.forEdges([&](uint32_t u, uint32_t v, float w) {
        ASSERT_TRUE(G2.hasEdge(u, v));
        EXPECT_EQ(w, G2.weight(u, v));
        EXPECT_EQ(G.edgeId(u, v), G2.edgeId(u, v));
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryUsesCompactV5TableWidth) {
    AdjListGraph<uint16_t, uint16_t> G(8, false, false);
    G.addEdge(0, 1);
    G.addEdge(2, 3);
    G.addEdge(4, 5);
    G.addEdge(6, 7);

    NetworkitBinaryWriter writer(8);
    const auto data = writer.writeToBuffer(G);

    uint64_t features;
    std::memcpy(&features, data.data() + 2 * sizeof(uint64_t), sizeof(uint64_t));
    const auto tableWidthCode = (features & nkbg::TABLE_WIDTH_MASK) >> nkbg::TABLE_WIDTH_SHIFT;
    EXPECT_EQ(nkbg::widthCode(1), tableWidthCode);
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryAdjListGraphVarintIntegerWeights) {
    using SmallGraph = AdjListGraph<uint16_t, int16_t>;

    SmallGraph G(4, true, false);
    G.addEdge(0, 1, -42);
    G.addEdge(2, 3, 327);

    NetworkitBinaryWriter writer(4);
    const auto data = writer.writeToBuffer(G);

    uint64_t features;
    std::memcpy(&features, data.data() + 2 * sizeof(uint64_t), sizeof(uint64_t));
    EXPECT_EQ(static_cast<uint64_t>(nkbg::WeightFormat::SIGNED_VARINT),
              (features & nkbg::WGHT_MASK) >> nkbg::WGHT_SHIFT);

    NetworkitBinaryReader reader;
    const SmallGraph G2 = reader.readFromBuffer<SmallGraph>(data);
    EXPECT_TRUE(G2.isWeighted());
    EXPECT_EQ(G.weight(0, 1), G2.weight(0, 1));
    EXPECT_EQ(G.weight(2, 3), G2.weight(2, 3));

    AdjListGraph<uint32_t, uint16_t> unsignedG(4, true, false);
    unsignedG.addEdge(0, 1, 42);
    unsignedG.addEdge(2, 3, 327);

    const auto unsignedData = writer.writeToBuffer(unsignedG);

    std::memcpy(&features, unsignedData.data() + 2 * sizeof(uint64_t), sizeof(uint64_t));
    EXPECT_EQ(static_cast<uint64_t>(nkbg::WeightFormat::VARINT),
              (features & nkbg::WGHT_MASK) >> nkbg::WGHT_SHIFT);

    const auto unsignedG2 = reader.readFromBuffer<AdjListGraph<uint32_t, uint16_t>>(unsignedData);
    EXPECT_TRUE(unsignedG2.isWeighted());
    EXPECT_EQ(unsignedG.weight(0, 1), unsignedG2.weight(0, 1));
    EXPECT_EQ(unsignedG.weight(2, 3), unsignedG2.weight(2, 3));
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryReadCompact) {
    NetworkitBinaryWriter writer(4);
    NetworkitBinaryReader reader;

    using CompactUnsignedGraph = AdjListGraph<uint16_t, uint8_t>;
    CompactUnsignedGraph unsignedG(300, true, false);
    unsignedG.addEdge(0, 299, 200);

    writer.write(unsignedG, "output/binary_compact_unsigned");
    const auto compactUnsigned = reader.readCompact("output/binary_compact_unsigned");
    ASSERT_TRUE(std::holds_alternative<CompactUnsignedGraph>(compactUnsigned));
    const auto &unsignedG2 = std::get<CompactUnsignedGraph>(compactUnsigned);
    EXPECT_EQ(unsignedG.numberOfNodes(), unsignedG2.numberOfNodes());
    EXPECT_EQ(unsignedG.weight(0, 299), unsignedG2.weight(0, 299));

    using SignedGraph = AdjListGraph<uint16_t, int16_t>;
    SignedGraph signedG(4, true, false);
    signedG.addEdge(0, 1, -42);
    signedG.addEdge(2, 3, 327);

    writer.write(signedG, "output/binary_compact_signed");
    const auto compactSigned = reader.readCompact("output/binary_compact_signed");
    using CompactSignedGraph = AdjListGraph<uint8_t, int16_t>;
    ASSERT_TRUE(std::holds_alternative<CompactSignedGraph>(compactSigned));
    const auto &signedG2 = std::get<CompactSignedGraph>(compactSigned);
    EXPECT_EQ(signedG.numberOfNodes(), signedG2.numberOfNodes());
    EXPECT_EQ(signedG.weight(0, 1), signedG2.weight(0, 1));
    EXPECT_EQ(signedG.weight(2, 3), signedG2.weight(2, 3));
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryTiny01Indexed) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/tiny_01.graph");
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);

    G.indexEdges();
    writer.write(G, "output/binary_tiny01");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_tiny01");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());

    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G2.edgeId(u, v), G.edgeId(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryKonect) {
    KONECTGraphReader reader2;
    Graph G = reader2.read("input/foodweb-baydry.konect");
    NetworkitBinaryWriter writer;

    writer.write(G, "output/binary_konect");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_konect");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), true);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G.weight(u, v), G2.weight(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryKonectInMemory) {
    KONECTGraphReader reader2;
    Graph G = reader2.read("input/foodweb-baydry.konect");
    NetworkitBinaryWriter writer;

    std::vector<uint8_t> data = writer.writeToBuffer(G);
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.readFromBuffer(data);
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), true);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G.weight(u, v), G2.weight(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryKonectIndexed) {
    KONECTGraphReader reader2;
    Graph G = reader2.read("input/foodweb-baydry.konect");
    G.indexEdges();
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);
    writer.write(G, "output/binary_konect");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_konect");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), true);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G.weight(u, v), G2.weight(u, v));
            ASSERT_EQ(G.edgeId(u, v), G2.edgeId(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryJazz) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/jazz.graph");

    NetworkitBinaryWriter writer;
    writer.write(G, "output/binary_jazz");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_jazz");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u) { G.forEdgesOf(u, [&](node v) { ASSERT_TRUE(G2.hasEdge(u, v)); }); });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryJazzIndexed) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/jazz.graph");
    G.indexEdges();

    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);
    writer.write(G, "output/binary_jazz");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_jazz");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u) { G.forEdgesOf(u, [&](node v) { ASSERT_TRUE(G2.hasEdge(u, v)); }); });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryWiki) {
    SNAPGraphReader reader2(true);
    Graph G = reader2.read("input/wiki-Vote.txt");
    NetworkitBinaryWriter writer;

    writer.write(G, "output/binary_wiki");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_wiki");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u) { G.forEdgesOf(u, [&](node v) { ASSERT_TRUE(G2.hasEdge(u, v)); }); });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryWikiIndexed) {
    SNAPGraphReader reader2(true);
    Graph G = reader2.read("input/wiki-Vote.txt");
    G.indexEdges();
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);

    writer.write(G, "output/binary_wiki");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_wiki");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u) { G.forEdgesOf(u, [&](node v) { ASSERT_TRUE(G2.hasEdge(u, v)); }); });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinarySignedWeights) {

    Graph G(10, true, false);
    int64_t weight = -1;
    for (count n = 0; n < G.numberOfNodes(); n++) {
        if (n != G.numberOfNodes() - 1)
            G.addEdge(n, n + 1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);
    writer.write(G, "output/binarySigned");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binarySigned");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G.weight(u, v), G2.weight(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinarySignedWeightsIndexed) {

    Graph G(10, true, false);
    G.indexEdges();
    int64_t weight = -1;
    for (count n = 0; n < G.numberOfNodes(); n++) {
        if (n != G.numberOfNodes() - 1)
            G.addEdge(n, n + 1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);
    writer.write(G, "output/binarySigned");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binarySigned");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G.weight(u, v), G2.weight(u, v));
            ASSERT_EQ(G.edgeId(u, v), G2.edgeId(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryFloatWeights) {

    Graph G(10, true, false);
    float weight = 987.654f;
    for (count n = 0; n < G.numberOfNodes(); n++) {
        if (n != G.numberOfNodes() - 1)
            G.addEdge(n, n + 1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);
    writer.write(G, "output/binaryFloats");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binaryFloats");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G.weight(u, v), G2.weight(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryFloatWeightsIndexed) {

    Graph G(10, true, false);
    G.indexEdges();
    float weight = 987.654f;
    for (count n = 0; n < G.numberOfNodes(); n++) {
        if (n != G.numberOfNodes() - 1)
            G.addEdge(n, n + 1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::AUTO_DETECT);
    writer.write(G, "output/binaryFloats");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binaryFloats");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u, v));
            ASSERT_EQ(G.weight(u, v), G2.weight(u, v));
        });
    });
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryUndirectedSelfLoops) {

    Graph G(5, false, false);
    G.addEdge(0, 0);
    G.addEdge(1, 1);
    G.addEdge(2, 2);
    G.addEdge(3, 3);
    G.addEdge(4, 4);
    NetworkitBinaryWriter writer;
    writer.write(G, "output/loops");
    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/loops");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G.numberOfSelfLoops(), G2.numberOfSelfLoops());
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryDirectedSelfLoops) {

    Graph G(5, false, true);
    G.addEdge(0, 0);
    G.addEdge(1, 1);
    G.addEdge(2, 2);
    G.addEdge(3, 3);
    G.addEdge(4, 4);
    NetworkitBinaryWriter writer;
    writer.write(G, "output/loops");
    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/loops");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G.numberOfSelfLoops(), G2.numberOfSelfLoops());
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryVarInt) {
    std::array<uint8_t, 10> buffer;

    // write defined values into buffer
    {
        uint8_t i = 0;
        for (auto &x : buffer)
            x = i++;
    }

    std::mt19937_64 gen{1};

    uint64_t checked_bits = 0;

    for (int bits = 0; bits < 64; ++bits) {
        auto min = uint64_t(1) << bits;
        auto max = 2 * min - 1;

        // special cases
        if (bits == 0) {
            min = 0;
            max = 0;
        } else if (bits == 64) {
            max = std::numeric_limits<uint64_t>::max();
        }

        std::uniform_int_distribution<uint64_t> distr{min, max};

        const auto nSamples = std::min<size_t>(10 * max + 2, 1000);

        for (size_t i = 0; i < nSamples; ++i) {
            // first two iterations test min/max values, all other random values
            const auto orig = [&] {
                if (i == 0)
                    return min;
                if (i == 1)
                    return max;
                return distr(gen);
            }();

            uint64_t valueRead;
            const auto bytesWritten = nkbg::varIntEncode(orig, buffer.data());
            const auto bytesRead = nkbg::varIntDecode(buffer.data(), valueRead);

            ASSERT_EQ(bytesWritten, bytesRead) << "bits=" << bits << ", i=" << i;
            ASSERT_EQ(valueRead, orig) << "bits=" << bits << ", i=" << i;
            ASSERT_GT(buffer[bytesWritten - 1], 0) << "bits=" << bits << ", i=" << i;

            for (size_t j = bytesWritten; j < buffer.size(); ++j)
                ASSERT_EQ(buffer[j], j);

            checked_bits |= orig;
        }
    }

    // make sure we touched each bit at least once
    ASSERT_EQ(checked_bits, std::numeric_limits<uint64_t>::max());
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryZigzag) {
    std::mt19937_64 gen(1);
    std::uniform_int_distribution<uint64_t> distr(0,
                                                  (std::numeric_limits<uint64_t>::max() >> 1) - 1);

    for (int i = 0; i < 10000; ++i) {
        auto check = [](int64_t value) {
            const auto encoded = nkbg::zigzagEncode(value);
            const auto decoded = nkbg::zigzagDecode(encoded);

            ASSERT_EQ(value, decoded);
            ASSERT_LE(encoded, 2u * static_cast<uint64_t>(std::abs(value)));
        };

        const auto x = distr(gen);
        check(x);
        check(-1 * x);
    }
}

TEST_F(NetworkitBinaryGTest, testNetworkitWriterNonContinuousNodesIds) {
    Graph G(20, true);
    G.removeNode(10);
    std::string path = "output/test.gt";
    NetworkitBinaryWriter{}.write(G, path);
    Graph GRead = NetworkitBinaryReader{}.read(path);
    EXPECT_EQ(G.numberOfNodes(), GRead.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), GRead.numberOfEdges());
    EXPECT_EQ(G.isDirected(), GRead.isDirected());
    EXPECT_EQ(G.isWeighted(), GRead.isWeighted());
    G.forNodes([&](node u) {
        G.forEdgesOf(u, [&](node v) {
            EXPECT_TRUE(GRead.hasEdge(u, v));
            EXPECT_DOUBLE_EQ(G.weight(u, v), GRead.weight(u, v));
        });
    });
}

MATCHER_P(GraphFeaturesEqual, expected, "Graph matches expected features") {
    return arg.numberOfNodes() == expected.numberOfNodes()
           && arg.numberOfEdges() == expected.numberOfEdges()
           && arg.isWeighted() == expected.isWeighted() && arg.isDirected() == expected.isDirected()
           && arg.hasEdgeIds() == expected.hasEdgeIds();
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryWriteReadEmptyGraph) {
    Graph graph(0, false, true);
    const std::filesystem::path path = "output/empty_graph.nkb";
    NetworkitBinaryWriter{}.write(graph, path.string());

    const Graph graph_read = NetworkitBinaryReader{}.read(path.string());
    EXPECT_THAT(graph_read, GraphFeaturesEqual(graph));
}

TEST_F(NetworkitBinaryGTest, testNetworkitBinaryWriteReadEmptyGraphWithIndexes) {
    Graph graph(0, true, false, true);
    const std::filesystem::path path = "output/empty_graph.nkb";
    NetworkitBinaryWriter{}.write(graph, path.string());

    const Graph graph_read = NetworkitBinaryReader{}.read(path.string());
    EXPECT_THAT(graph_read, GraphFeaturesEqual(graph));
}

TEST_F(NetworkitBinaryGTest, testWriteReadNonContinuousUndirectedWeighted) {
    Graph graph(8, true, false);

    graph.removeNode(1);
    graph.removeNode(4);
    graph.removeNode(6);

    graph.addEdge(0, 2, 2.0);
    graph.addEdge(0, 3, 1.5);
    graph.addEdge(2, 3, 3.25);
    graph.addEdge(2, 5, 0.75);
    graph.addEdge(3, 7, 4.0);
    graph.addEdge(5, 7, 2.5);

    // Edge that skips removed IDs
    graph.addEdge(0, 7, 9.0);

    std::string path = "output/writereadnoncont_weighted_undirected.nkb";

    NetworkitBinaryWriter writer;
    writer.write(graph, path);

    NetworkitBinaryReader reader;
    Graph readGraph = reader.read(path);

    EXPECT_EQ(readGraph.numberOfNodes(), graph.numberOfNodes());
    EXPECT_EQ(readGraph.upperNodeIdBound(), graph.upperNodeIdBound());
    EXPECT_EQ(readGraph.isDirected(), graph.isDirected());
    EXPECT_EQ(readGraph.isWeighted(), graph.isWeighted());
    EXPECT_EQ(readGraph.numberOfEdges(), graph.numberOfEdges());

    for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
        EXPECT_EQ(readGraph.hasNode(u), graph.hasNode(u)) << "u=" << u;
    }

    for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
        if (!graph.hasNode(u))
            continue;

        EXPECT_EQ(readGraph.degree(u), graph.degree(u)) << "deg u=" << u;

        std::vector<node> nbrG, nbrR;
        graph.forNeighborsOf(u, [&](node v) { nbrG.push_back(v); });
        readGraph.forNeighborsOf(u, [&](node v) { nbrR.push_back(v); });

        std::ranges::sort(nbrG);
        std::ranges::sort(nbrR);
        EXPECT_EQ(nbrR, nbrG) << "neighbors u=" << u;
    }

    // Adjacency equality per node
    for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
        std::vector<node> nbrG, nbrR;

        graph.forNeighborsOf(u, [&](node v) { nbrG.push_back(v); });
        readGraph.forNeighborsOf(u, [&](node v) { nbrR.push_back(v); });

        std::ranges::sort(nbrG);
        std::ranges::sort(nbrR);

        EXPECT_EQ(nbrR, nbrG) << "neighbors u=" << u;
    }

    const std::vector<WeightedEdge> edges = {
        {0, 2, 2.0}, {0, 3, 1.5}, {2, 3, 3.25}, {2, 5, 0.75}, {3, 7, 4.0}, {5, 7, 2.5}, {0, 7, 9.0},
    };

    for (const auto &e : edges) {
        EXPECT_TRUE(readGraph.hasEdge(e.u, e.v)) << "edge " << e.u << "-" << e.v;
        EXPECT_TRUE(readGraph.hasEdge(e.v, e.u)) << "edge " << e.v << "-" << e.u;

        EXPECT_DOUBLE_EQ(readGraph.weight(e.u, e.v), e.weight) << "w(" << e.u << "," << e.v << ")";
        EXPECT_DOUBLE_EQ(readGraph.weight(e.v, e.u), e.weight) << "w(" << e.v << "," << e.u << ")";
    }

    // Sanity check: removed nodes have no incident edges
    for (node u : {1u, 4u, 6u}) {
        EXPECT_FALSE(readGraph.hasNode(u));
    }
}

TEST_F(NetworkitBinaryGTest, testWriteReadNonContinuousDirectedEdgesIndexed) {
    Graph graph(10, false, true);

    graph.removeNode(1);
    graph.removeNode(4);
    graph.removeNode(8);

    graph.addEdge(0, 2);
    graph.addEdge(2, 0);
    graph.addEdge(0, 3);
    graph.addEdge(3, 5);
    graph.addEdge(5, 7);
    graph.addEdge(7, 9);
    graph.addEdge(9, 2);
    graph.addEdge(6, 7);
    graph.addEdge(6, 9);
    graph.addEdge(2, 7);
    graph.addEdge(9, 0);

    graph.indexEdges();

    std::string path = "output/writereadnoncont_directed_edgeids.nkb";

    NetworkitBinaryWriter writer;
    writer.write(graph, path);

    NetworkitBinaryReader reader;
    Graph readGraph = reader.read(path);

    EXPECT_EQ(readGraph.numberOfNodes(), graph.numberOfNodes());
    EXPECT_EQ(readGraph.upperNodeIdBound(), graph.upperNodeIdBound());
    EXPECT_EQ(readGraph.isDirected(), graph.isDirected());
    EXPECT_EQ(readGraph.isWeighted(), graph.isWeighted());
    EXPECT_EQ(readGraph.numberOfEdges(), graph.numberOfEdges());

    for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
        EXPECT_EQ(readGraph.hasNode(u), graph.hasNode(u)) << "u=" << u;
    }

    // Adjacency equality (strong)
    for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
        std::vector<node> nbrG, nbrR;
        graph.forNeighborsOf(u, [&](node v) { nbrG.push_back(v); });
        readGraph.forNeighborsOf(u, [&](node v) { nbrR.push_back(v); });

        std::ranges::sort(nbrG);
        std::ranges::sort(nbrR);
        EXPECT_EQ(nbrR, nbrG) << "neighbors u=" << u;
    }

    // Edge ID checks
    EXPECT_TRUE(graph.hasEdgeIds());
    EXPECT_TRUE(readGraph.hasEdgeIds());

    graph.forEdges([&](node u, node v) {
        const auto idG = graph.edgeId(u, v);
        const auto idR = readGraph.edgeId(u, v);
        EXPECT_EQ(idR, idG) << "edgeId(" << u << "->" << v << ")";
    });

    for (edgeid id = 0; id < graph.numberOfEdges(); ++id) {
        auto [u, v] = readGraph.edgeById(id);

        EXPECT_TRUE(readGraph.hasEdge(u, v)) << "edgeById(" << id << ") gives non-edge";
        EXPECT_EQ(readGraph.edgeId(u, v), id) << "edgeId(edgeById(" << id << "))";
    }
}

} /* namespace NetworKit */
