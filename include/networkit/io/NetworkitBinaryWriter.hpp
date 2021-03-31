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

enum class NetworkitBinaryWeights {
    none,
    unsignedFormat,
    signedFormat,
    doubleFormat,
    floatFormat,
    autoDetect
};

enum class NetworkitBinaryEdgeIDs { noEdgeIDs, writeEdgeIDs, autoDetect };

/**
 * @ingroup io
 *
 * Writes a graph in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */
class NetworkitBinaryWriter final : public GraphWriter {

public:
    NetworkitBinaryWriter(uint64_t chunks = 32,
                          NetworkitBinaryWeights weightsType = NetworkitBinaryWeights::autoDetect,
                          NetworkitBinaryEdgeIDs edgeIndex = NetworkitBinaryEdgeIDs::autoDetect);

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
