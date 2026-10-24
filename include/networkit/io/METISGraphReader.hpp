/*
 * METISGraphReader.hpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_GRAPH_READER_HPP_
#define NETWORKIT_IO_METIS_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the METIS file format documented in [1]
 *
 * The file does not store the node and edge weight types: they are given by @a GraphT.
 * Reading fails if the number of nodes does not fit into GraphT::NodeT, or if GraphT::EdgeWeightT
 * is integral and an edge weight in the file is not an integer representable by it.
 *
 * [1] http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf
 */
template <class GraphT>
class GenericMETISGraphReader final : public GraphReaderBase<GraphT> {
public:
    GenericMETISGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph file.
     *
     * @param[in]  path  file path
     *
     * @param[out]  the graph read from file
     */
    GraphT read(std::string_view path) override;
};

using METISGraphReader = GenericMETISGraphReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/METISGraphReaderImpl.hpp>

#endif // NETWORKIT_IO_METIS_GRAPH_READER_HPP_
