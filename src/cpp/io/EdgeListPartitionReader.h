#ifndef EDGELISTPARTITION_H_
#define EDGELISTPARTITION_H_

#include <fstream>

#include "../structures/Partition.h"
#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"


namespace NetworKit {

class EdgeListPartitionReader {

public:

	EdgeListPartitionReader(node firstNode=1);

	virtual ~EdgeListPartitionReader();

	/**
	 * Read a clustering from a file. File format:
	 * 		A list of the nodes and their membership (memberships are labelled by integer numbers >=1).
	 *
	 * @param[in]	path	path to file
	 */
	virtual Partition read(std::string path);

	/**
	* @param[in]	firstNode	index of the first node in the file
	*/
	node firstNode;
};

} /* namespace NetworKit */
#endif /* PARTITIONREADER_H_ */
