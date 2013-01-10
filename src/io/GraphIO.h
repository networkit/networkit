/*
 * GraphIO.h
 *
 *  Created on: 09.01.2013
 *      Author: cls
 */

#ifndef GRAPHIO_H_
#define GRAPHIO_H_

#include <string>
#include <iostream>
#include <fstream>

#include "../graph/Graph.h"
#include "../aux/log.h"

namespace EnsembleClustering {

class GraphIO {
public:

	GraphIO();

	virtual ~GraphIO();

	virtual void toEdgeList(Graph& G, std::string path);
};

} /* namespace EnsembleClustering */
#endif /* GRAPHIO_H_ */
