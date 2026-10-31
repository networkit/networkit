/*
 * GraphReader.hpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_GRAPH_READER_HPP_
#define NETWORKIT_IO_GRAPH_READER_HPP_

#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Type-independent part of all graph readers. Shared by all instantiations of
 * GraphReaderBase, so that e.g. GraphReader::DISCARD_EDGES has a single type.
 */
class GraphReaderOptions {
public:
    // How to deal with multi-edges
    enum MultipleEdgesHandling {
        DISCARD_EDGES,  // Reads and selects the first edge which occurs and discards all following
        SUM_WEIGHTS_UP, // If an edge occurs again, the weight of it is added to the existing edge
        KEEP_MINIMUM_WEIGHT // The edge with the lowest weight is kept
    };
};

/**
 * @ingroup io
 * Abstract base class for graph readers.
 */
template <class GraphT>
class GraphReaderBase : public GraphReaderOptions {
public:
    virtual ~GraphReaderBase() = default;

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    virtual GraphT read(std::string_view path) = 0;
};

using GraphReader = GraphReaderBase<Graph>;

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GRAPH_READER_HPP_
