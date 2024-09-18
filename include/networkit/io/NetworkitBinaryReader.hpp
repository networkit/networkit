/*
 * NetworkitBinaryReader.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_READER_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_READER_HPP_

#include <cstring>
#include <string>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>
#include <networkit/io/MemoryMappedFile.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Reads a graph written in the custom Networkit binary format.
 * Note that there are multiple versions of the Networkit binary format.
 * This Reader can read files that are written with format version 2 and 3.
 * Format version 3 was released with Networkit 9.1 (December 2021).
 * Further information can be found here:
 * https://github.com/networkit/networkit/blob/master/networkit/cpp/io/NetworkitBinaryGraph.md
 */

class NetworkitBinaryReader final : public GraphReader {

public:
    NetworkitBinaryReader(){};

    Graph read(std::string_view path) override;
    Graph readFromBuffer(const std::vector<uint8_t> &data);

private:
    count nodes;
    count chunks;
    bool directed;
    bool weighted;
    bool indexed;
    count version;

    template <class T>
    Graph readData(const T &source);
};
} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_READER_HPP_
