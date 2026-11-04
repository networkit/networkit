/*
 * EdgeListWriter.hpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_EDGE_LIST_WRITER_HPP_
#define NETWORKIT_IO_EDGE_LIST_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * A writer for the edge list format. The output will contain one edge per line,
 * in the format fromNodeSEPARATORtoNode, where separator can be specified by
 * the user.
 */
template <class GraphT>
class GenericEdgeListWriter final : public GraphWriterBase<GraphT> {
    using NodeT = typename GraphT::NodeT;

public:
    GenericEdgeListWriter() = default; // nullary constructor for Python shell

    /**
     * @param[in]  separator  character used to separate nodes in an edge line
     * @param[in]  firstNode  index of the first node in the file
     * @param[in]  bothDirections  for undirected graphs: if every edge shall be written in both
     * directions (default: false)
     */
    GenericEdgeListWriter(char separator, NodeT firstNode, bool bothDirections = false);

    /**
     * Write the graph to a file.
     * @param[in]  G    the graph
     * @param[in]  path  the output file path
     */
    void write(const GraphT &G, std::string_view path) override;

private:
    char separator; //!< character separating nodes in an edge line
    NodeT firstNode;
    bool bothDirections;
};

using EdgeListWriter = GenericEdgeListWriter<Graph>;

} /* namespace NetworKit */

#include <networkit/io/EdgeListWriterImpl.hpp>

#endif // NETWORKIT_IO_EDGE_LIST_WRITER_HPP_
