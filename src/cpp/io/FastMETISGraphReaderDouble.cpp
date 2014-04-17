/*
 * FastMETISGraphReader.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "FastMETISGraphReaderDouble.h"

namespace NetworKit {

FastMETISGraphReaderDouble::FastMETISGraphReaderDouble() {
	// TODO Auto-generated constructor stub

}

FastMETISGraphReaderDouble::~FastMETISGraphReaderDouble() {
	// TODO Auto-generated destructor stub
}

Graph FastMETISGraphReaderDouble::read(std::string path) {
	FastMETISParserDouble parser;
	return parser.parse(path);
}

} /* namespace NetworKit */
