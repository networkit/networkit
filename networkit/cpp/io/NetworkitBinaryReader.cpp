/*
 * NetworkitBinaryReader.cpp
 *
 * Author: Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 */

#include <networkit/io/NetworkitBinaryReader.hpp>

namespace NetworKit {

Graph NetworkitBinaryReader::read(std::string_view path) {
    return read<Graph>(path);
}

Graph NetworkitBinaryReader::readFromBuffer(const std::vector<uint8_t> &data) {
    return readFromBuffer<Graph>(data);
}

} // namespace NetworKit
