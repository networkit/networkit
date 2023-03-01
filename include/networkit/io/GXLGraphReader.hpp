#ifndef NETWORKIT_IO_GXL_GRAPH_READER_HPP_
#define NETWORKIT_IO_GXL_GRAPH_READER_HPP_

#include "networkit/io/GraphReader.hpp"

namespace NetworKit {

/**
 * @ingroup io
 */
class GXLGraphReader final : public GraphReader {
public:
    GXLGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph file.
     *
     * @param[in]  path  file path
     *
     * @param[out]  the graph read from file
     */
    Graph read(const std::string &path) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GXL_GRAPH_READER_HPP_