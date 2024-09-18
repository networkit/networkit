#ifndef NETWORKIT_IO_MTX_GRAPH_READER_HPP_
#define NETWORKIT_IO_MTX_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Reader for the matrix market file format documented in
 * https://networkrepository.com/mtx-matrix-market-format.html
 *
 * Does not allow complex fields.
 *
 */
class MTXGraphReader final : public GraphReader {
public:
    MTXGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph.
     *
     * @param path
     * @return Graph
     */
    Graph read(std::string_view path) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_MTX_GRAPH_READER_HPP_
