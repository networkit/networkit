#ifndef EDGELISTCLUSTERING_H_
#define EDGELISTCLUSTERING_H_

#include <fstream>

#include "../structures/Partition.h"
#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"


namespace NetworKit {

class EdgeListClusteringReader {

public:

	EdgeListClusteringReader(node firstNode=1);

	virtual ~EdgeListClusteringReader();

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
#endif /* CLUSTERINGREADER_H_ */
