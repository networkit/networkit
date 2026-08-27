/*
 * NetworkitBinaryWriter.cpp
 *
 * @author Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

#include <networkit/io/NetworkitBinaryWriter.hpp>

namespace NetworKit {

NetworkitBinaryWriter::NetworkitBinaryWriter(uint64_t chunks, NetworkitBinaryWeights weightsType,
                                             NetworkitBinaryEdgeIDs edgeIndex)
    : chunks(chunks), weightsType(weightsType), edgeIndex(edgeIndex) {}

void NetworkitBinaryWriter::write(const Graph &G, std::string_view path) {
    write<Graph>(G, path);
}

std::vector<uint8_t> NetworkitBinaryWriter::writeToBuffer(const Graph &G) {
    return writeToBuffer<Graph>(G);
}

} // namespace NetworKit
