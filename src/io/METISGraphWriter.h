/*
 * METISGraphWriter.h
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef METISGRAPHWRITER_H_
#define METISGRAPHWRITER_H_

#include <fstream>

#include "GraphWriter.h"

namespace EnsembleClustering {

class METISGraphWriter: public EnsembleClustering::GraphWriter {

public:

	METISGraphWriter();

	virtual ~METISGraphWriter();

	virtual void write(Graph& G, std::string path);
};

} /* namespace EnsembleClustering */
#endif /* METISGRAPHWRITER_H_ */
