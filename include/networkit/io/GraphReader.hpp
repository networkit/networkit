/*
 * GraphReader.h
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_IO_GRAPH_READER_HPP_
#define NETWORKIT_IO_GRAPH_READER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/StringTools.hpp>


namespace NetworKit {

/*
* Enum for defining handling of graph formats with multiple same edges.
*/



/**
 * @ingroup io
 * Abstract base class for graph readers.
 */
class GraphReader {
public:

    enum MultipleEdgesHandling
    {
         DISCARD_EDGES, //Reads and selects the first edge which occurs and discards all following
         SUM_WEIGHTS_UP, //If an edge occurs again, the weight of it is added to the existing edge
         KEEP_MINIUM_WEIGHT //The edge with the lowest weight is kept
    };

    virtual ~GraphReader() = default;

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]	path	input file path
     */
    virtual Graph read(const std::string& path) = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GRAPH_READER_HPP_
