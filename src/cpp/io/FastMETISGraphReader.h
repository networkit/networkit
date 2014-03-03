/*
 * FastMETISGraphReader.h
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef FASTMETISGRAPHREADER_H_
#define FASTMETISGRAPHREADER_H_

#include <fstream>

#include "GraphReader.h"
#include "FastMETISParser.h"

namespace NetworKit {

class FastMETISGraphReader: public NetworKit::GraphReader {
public:
	FastMETISGraphReader();
	virtual ~FastMETISGraphReader();
	virtual Graph read(std::string path);
};

} /* namespace NetworKit */
#endif /* FASTMETISGRAPHREADER_H_ */
