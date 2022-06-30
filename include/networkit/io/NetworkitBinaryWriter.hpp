/*
 * NetworkitBinaryWriter.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_WRITER_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_WRITER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

enum class NetworkitBinaryWeights : int {
    NONE,
    UNSIGNED_FORMAT,
    SIGNED_FORMAT,
    DOUBLE_FORMAT,
    FLOAT_FORMAT,
    AUTO_DETECT,
    none = NONE, // this + following added for backwards compatibility
    unsignedFormat = UNSIGNED_FORMAT,
    signedFormat = SIGNED_FORMAT,
    doubleFormat = DOUBLE_FORMAT,
    floatFormat = FLOAT_FORMAT,
    autoDetect = AUTO_DETECT
};

enum class NetworkitBinaryEdgeIDs : int {
    NO_EDGE_IDS,
    WIRTE_EDGE_IDS,
    AUTO_DETECT,
    noEdgeIDs = NO_EDGE_IDS, // this + following added for backwards compatibility
    writeEdgeIDs = WIRTE_EDGE_IDS,
    autoDetect = AUTO_DETECT
};

/**
 * @ingroup io
 *
 * Writes a graph in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */
class NetworkitBinaryWriter final : public GraphWriter {

public:
    NetworkitBinaryWriter(uint64_t chunks = 32,
                          NetworkitBinaryWeights weightsType = NetworkitBinaryWeights::AUTO_DETECT,
                          NetworkitBinaryEdgeIDs edgeIndex = NetworkitBinaryEdgeIDs::AUTO_DETECT);

    void write(const Graph &G, const std::string &path) override;
    std::vector<uint8_t> writeToBuffer(const Graph &G);

private:
    count chunks;
    NetworkitBinaryWeights weightsType;
    NetworkitBinaryEdgeIDs edgeIndex;
    bool preserveEdgeIndex;

    template <class T>
    void writeData(T &outStream, const Graph &G);
};

} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_WRITER_HPP_
