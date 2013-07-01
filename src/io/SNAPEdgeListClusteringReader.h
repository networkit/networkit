/*
 * SNAPEdgeListClusteringReader.h
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */

#ifndef SNAPEDGELISTCLUSTERINGREADER_H_
#define SNAPEDGELISTCLUSTERINGREADER_H_

#include <unordered_set>
#include <vector>
#include <fstream>

#include "../clustering/Clustering.h"
#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

class SNAPEdgeListClusteringReader {
public:
	SNAPEdgeListClusteringReader();
	virtual ~SNAPEdgeListClusteringReader();

	virtual std::vector<std::set<node>>  read(std::string path);

};

} /* namespace NetworKit */
#endif /* SNAPEDGELISTCLUSTERINGREADER_H_ */
