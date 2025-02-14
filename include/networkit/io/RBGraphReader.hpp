/*
 * RBGraphReader.hpp
 *
 *  Created on: 16.10.2024
 *      Author: bernlu
 */

#ifndef NETWORKIT_IO_RB_GRAPH_READER_HPP_
#define NETWORKIT_IO_RB_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the Rutherford Boeing (RB) matrix file format as described in
 * http://sparse-files.engr.tamu.edu/files/DOC/rb.pdf.
 *
 * @note currently the reader only supports compressed column format for real, integer, or pattern
 * data types.
 */
class RBGraphReader final : public GraphReader {
public:
    RBGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph.
     *
     * @param path
     * @return Graph
     */
    Graph read(std::string_view path) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_RB_GRAPH_READER_HPP_
