/*
 * SNAPEdgeListPartitionReader.hpp
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */

#ifndef NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_
#define NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_

#include <fstream>
#include <set>
#include <unordered_map>
#include <vector>

#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reads the clustering files from the SNAP collection.
 */
class SNAPEdgeListPartitionReader final {
public:
    Cover read(std::string path, std::unordered_map<node,node>& mapNodeIds, Graph& G);

};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_
