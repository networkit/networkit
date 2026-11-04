/*
 * GraphWriter.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_GRAPH_WRITER_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Abstract base class for graph writers.
 */
template <class GraphT>
class GraphWriterBase {
public:
    virtual ~GraphWriterBase() = default;

    virtual void write(const GraphT &G, std::string_view path) = 0;
};

using GraphWriter = GraphWriterBase<Graph>;

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GRAPH_WRITER_HPP_
