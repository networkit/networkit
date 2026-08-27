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
 * This Reader can read files that are written with format versions 2 through 5.
 * Format version 3 was released with Networkit 9.1 (December 2021).
 * Further information can be found here:
 * https://github.com/networkit/networkit/blob/master/networkit/cpp/io/NetworkitBinaryGraph.md
 *
 * 	Note
 *  ----
 * Version 5 supports reading templated AdjListGraph instances and compact chunk metadata tables.
 *
 */

class NetworkitBinaryReader final : public GraphReader {

public:
    NetworkitBinaryReader() {};

    Graph read(std::string_view path) override;
    template <class GraphT>
    GraphT read(std::string_view path);
    Graph readFromBuffer(const std::vector<uint8_t> &data);
    template <class GraphT>
    GraphT readFromBuffer(const std::vector<uint8_t> &data);

private:
    count nodes;
    count chunks;
    bool directed;
    bool weighted;
    bool indexed;
    count version;
    uint8_t tableWidth;

    template <class GraphT, class T>
    GraphT readData(const T &source);
};
} // namespace NetworKit

#include <networkit/io/NetworkitBinaryReaderImpl.hpp>

#endif // NETWORKIT_IO_NETWORKIT_BINARY_READER_HPP_
