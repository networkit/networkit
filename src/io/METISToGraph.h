/*
 * METISToGraph.h
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#ifndef METISTOGRAPH_H_
#define METISTOGRAPH_H_

#include <string>

#include "../graph/Graph.h"

namespace EnsembleClustering {


/**
 * This class provides a user interface for reading a METIS graph file and returning a
 * Graph object.
 */
class METISToGraph {

public:

	METISToGraph();

	virtual ~METISToGraph();

	virtual Graph read(std::string graphPath);
};

} /* namespace EnsembleClustering */
#endif /* METISTOSTINGER_H_ */
