/*
 * GraphWriter.h
 *
 *  Created on: 30.01.2013
 *      Author: cls
 */

#ifndef GRAPHWRITER_H_
#define GRAPHWRITER_H_

#include "../graph/Graph.h"

namespace EnsembleClustering {

class GraphWriter {
public:

	GraphWriter();

	virtual ~GraphWriter();

	virtual void write(Graph& G, std::string path) = 0;
};

} /* namespace EnsembleClustering */
#endif /* GRAPHWRITER_H_ */
