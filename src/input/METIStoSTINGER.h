/*
 * METIStoSTINGER.h
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#ifndef METISTOSTINGER_H_
#define METISTOSTINGER_H_

#include <string>

#include "../graph/Graph.h"

namespace EnsembleClustering {


/**
 * This class provides a user interface for reading a METIS graph file and returning a
 * STINGER-based graph object.
 */
class METIStoSTINGER {

public:

	METIStoSTINGER();

	virtual ~METIStoSTINGER();

	virtual Graph* read(std::string graphPath);
};

} /* namespace EnsembleClustering */
#endif /* METISTOSTINGER_H_ */
