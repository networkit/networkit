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

namespace NetworKit {

/**
 * @ingroup io
 */
class METISGraphWriter: public NetworKit::GraphWriter {

public:

	METISGraphWriter() = default;

	virtual void write(Graph& G, const std::string& path) override;

	virtual void write(Graph& G, bool weighted, std::string path);
};

} /* namespace NetworKit */
#endif /* METISGRAPHWRITER_H_ */
