/*
 * FastMETISGraphReader.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "FastMETISGraphReader.h"

namespace NetworKit {

FastMETISGraphReader::FastMETISGraphReader() {
	// TODO Auto-generated constructor stub

}

FastMETISGraphReader::~FastMETISGraphReader() {
	// TODO Auto-generated destructor stub
}

Graph FastMETISGraphReader::read(std::string path) {
	FastMETISParser parser;
	return parser.parse(path);
}

} /* namespace NetworKit */
