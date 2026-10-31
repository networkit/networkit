/*
 * METISGraphWriter.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_METIS_GRAPH_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
template <class GraphT>
class GenericMETISGraphWriter final : public GraphWriterBase<GraphT> {

public:
    GenericMETISGraphWriter() = default;

    void write(const GraphT &G, std::string_view path) override;

    void write(const GraphT &G, bool weighted, std::string_view path);
};

using METISGraphWriter = GenericMETISGraphWriter<Graph>;

} /* namespace NetworKit */

#include <networkit/io/METISGraphWriterImpl.hpp>

#endif // NETWORKIT_IO_METIS_GRAPH_WRITER_HPP_
