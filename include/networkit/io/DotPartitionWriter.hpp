/*
 * DotPartitionWriter.h
 */

#ifndef NETWORKIT_IO_DOT_PARTITION_WRITER_HPP_
#define NETWORKIT_IO_DOT_PARTITION_WRITER_HPP_

#include <map>

#include <networkit/structures/Partition.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class DotPartitionWriter {
public:
    virtual void write(Graph& graph, Partition& zeta, std::string path) const;

    virtual std::map<index, double> createHueMap(Graph &graph, Partition& zeta) const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_DOT_PARTITION_WRITER_HPP_
