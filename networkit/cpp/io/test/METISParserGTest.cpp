#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISParser.hpp>

#include "IOTestUtils.hpp"

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::HasSubstr;
using ::testing::Pair;
using ::testing::ThrowsMessage;
using TestUtils::TemporaryDirectory;

// A METIS file with the given header and adjacency line in a temporary directory, and a parser
// of it that has read the header. The directory outlives the parser, since open files cannot be
// removed on all platforms.
struct TemporaryMETISFile {
    TemporaryMETISFile(const std::string &header, const std::string &line)
        : parser(directory.writeFile("graph.metis", header + "\n" + line + "\n")) {
        parser.getHeader();
    }

    const TemporaryDirectory directory;
    METISParser parser;
};

// Parses the weighted adjacency line of the first node of a METIS file with the given lines.
template <class NodeT, class EdgeWeightT>
std::vector<std::pair<NodeT, EdgeWeightT>> parseFirstWeightedLine(const std::string &header,
                                                                  const std::string &line) {
    return TemporaryMETISFile(header, line).parser.getNextWithWeights<NodeT, EdgeWeightT>();
}

// Parses the unweighted adjacency line of the first node of a METIS file with the given lines.
template <class NodeT>
std::vector<NodeT> parseFirstLine(const std::string &header, const std::string &line) {
    return TemporaryMETISFile(header, line).parser.getNext<NodeT>();
}

TEST(METISParserGTest, testParsesIntegerWeightsExactly) {
    // 2^53 + 1 is not representable by double, so a detour through double would round it.
    constexpr uint64_t notRepresentableByDouble = (uint64_t{1} << 53) + 1;
    const auto unsignedAdjacencies = parseFirstWeightedLine<node, uint64_t>(
        "2 1 1", "2 " + std::to_string(notRepresentableByDouble));
    EXPECT_THAT(unsignedAdjacencies, ElementsAre(Pair(node{2}, notRepresentableByDouble)));
    constexpr int64_t minInt64 = std::numeric_limits<int64_t>::min();
    const auto signedAdjacencies =
        parseFirstWeightedLine<node, int64_t>("2 1 1", "2 " + std::to_string(minInt64));
    EXPECT_THAT(signedAdjacencies, ElementsAre(Pair(node{2}, minInt64)));
}

TEST(METISParserGTest, testParsesFloatingPointWeights) {
    const auto floatAdjacencies = parseFirstWeightedLine<node, float>("3 2 1", "  2 0.5\t3 1.25  ");
    EXPECT_THAT(floatAdjacencies, ElementsAre(Pair(node{2}, 0.5f), Pair(node{3}, 1.25f)));
    const auto doubleAdjacencies = parseFirstWeightedLine<node, double>("3 2 1", "2 0.5 3 1.25");
    EXPECT_THAT(doubleAdjacencies, ElementsAre(Pair(node{2}, 0.5), Pair(node{3}, 1.25)));
}

TEST(METISParserGTest, testThrowsOnFractionalIntegerWeight) {
    const auto parse = []() -> void { parseFirstWeightedLine<node, int32_t>("2 1 1", "2 1.5"); };
    EXPECT_THAT(parse, ThrowsMessage<std::runtime_error>(
                           HasSubstr("METIS value 1.5 is not valid for the requested type")));
}

TEST(METISParserGTest, testThrowsOnNegativeUnsignedWeight) {
    const auto parse = []() -> void { parseFirstWeightedLine<node, uint32_t>("2 1 1", "2 -1"); };
    EXPECT_THAT(parse, ThrowsMessage<std::runtime_error>(
                           HasSubstr("METIS value -1 is not valid for the requested type")));
}

TEST(METISParserGTest, testThrowsOnOutOfRangeWeight) {
    const auto parse = []() -> void { parseFirstWeightedLine<node, uint8_t>("2 1 1", "2 256"); };
    EXPECT_THAT(parse, ThrowsMessage<std::runtime_error>(
                           HasSubstr("METIS value 256 is out of range for the requested type")));
}

TEST(METISParserGTest, testThrowsOnMissingWeight) {
    const auto parse = []() -> void { parseFirstWeightedLine<node, double>("3 2 1", "2 0.5 3"); };
    EXPECT_THAT(parse, ThrowsMessage<std::runtime_error>(
                           HasSubstr("METIS line ended where a value was expected")));
}

TEST(METISParserGTest, testThrowsOnOutOfRangeNodeId) {
    // 2^64 does not fit into node, so it can only be spelled as a string.
    constexpr const char *tooLargeNodeId = "18446744073709551616";
    EXPECT_THAT(
        []() -> void { parseFirstLine<node>("2 1", tooLargeNodeId); },
        ThrowsMessage<std::runtime_error>(HasSubstr(std::string("METIS value ") + tooLargeNodeId
                                                    + " is out of range for the requested type")));
}

TEST(METISParserGTest, testParsesNodeIdsIntoRequestedType) {
    EXPECT_THAT(parseFirstLine<int32_t>("3 2", "2 3"), ElementsAre(int32_t{2}, int32_t{3}));
    const auto weightedAdjacencies =
        parseFirstWeightedLine<int32_t, double>("3 2 1", "2 0.5 3 1.25");
    EXPECT_THAT(weightedAdjacencies, ElementsAre(Pair(int32_t{2}, 0.5), Pair(int32_t{3}, 1.25)));
}

TEST(METISParserGTest, testThrowsOnNodeIdOutOfRangeOfRequestedType) {
    // 2^31 does not fit into int32_t.
    EXPECT_THAT([]() -> void { parseFirstLine<int32_t>("2 1", "2147483648"); },
                ThrowsMessage<std::runtime_error>(
                    HasSubstr("METIS value 2147483648 is out of range for the requested type")));
    const auto parse = []() -> void {
        parseFirstWeightedLine<int32_t, double>("2 1 1", "2147483648 1");
    };
    EXPECT_THAT(parse, ThrowsMessage<std::runtime_error>(HasSubstr(
                           "METIS value 2147483648 is out of range for the requested type")));
}

TEST(METISParserGTest, testReaderThrowsOnMalformedWeightedLine) {
    // Malformed weighted lines used to be skipped with an error message, dropping their edges.
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("graph.metis", "2 1 1\n2 x\n1 x\n");
    EXPECT_THAT([&]() -> void { METISGraphReader{}.read(path); },
                ThrowsMessage<std::runtime_error>(
                    HasSubstr("METIS value x is not valid for the requested type")));
}

} // namespace
} // namespace NetworKit
