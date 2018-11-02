#ifndef EDGELISTPARTITION_H_
#define EDGELISTPARTITION_H_

#include <fstream>

#include "../structures/Partition.h"
#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"


namespace NetworKit {

/**
 * @ingroup io
 */
class EdgeListPartitionReader {

public:

	/**
	 * Constructs the EdgeListPartitionReader class with @a firstNode as the index of the first node in the file.
	 * @param[in]	firstNode	Index of the first node in the file.
	 * @param[in]	sepChar		The separator between two elements
	 */
	EdgeListPartitionReader(node firstNode=1, char sepChar='\t');

	/**
	 * Read a clustering from a file. File format:
	 * 		A list of the nodes and their membership (memberships are labelled by integer numbers >=1).
	 *
	 * @param[in]	path	Path to file.
	 * @return The clustering contained in the file at @a path.
	 */
	virtual Partition read(std::string path);


	node firstNode;
	char sepChar;
};

} /* namespace NetworKit */
#endif /* PARTITIONREADER_H_ */
