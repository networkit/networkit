/*
 * FastMETISGraphReader.h
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef FASTMETISGRAPHREADERDOUBLE_H_
#define FASTMETISGRAPHREADERDOUBLE_H_

#include <fstream>

#include "GraphReader.h"
#include "FastMETISParserDouble.h"

namespace NetworKit {

class FastMETISGraphReaderDouble: public NetworKit::GraphReader {
public:
	FastMETISGraphReaderDouble();
	virtual ~FastMETISGraphReaderDouble();
	virtual Graph read(std::string path);
};

} /* namespace NetworKit */
#endif /* FASTMETISGRAPHREADER_H_ */
