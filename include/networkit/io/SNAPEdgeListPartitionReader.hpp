/*
 * SNAPEdgeListPartitionReader.h
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */

#ifndef NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_
#define NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_

#include <set>
#include <vector>
#include <fstream>
#include <unordered_map>

#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reads the clustering files from the SNAP collection.
 */
class SNAPEdgeListPartitionReader {
public:
    virtual Cover read(std::string path, std::unordered_map<node,node>& mapNodeIds, Graph& G);

//	virtual Partition readWithInfo(std::string path, count nNodes);

};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_
