/*
 * DGSReader.cpp
 *
 *  Created on: 01.06.2013
 *      Author: cls
 */

#include "DGSReader.h"

namespace NetworKit {

DGSReader::DGSReader() {
	// TODO Auto-generated constructor stub

}

DGSReader::~DGSReader() {
	// TODO Auto-generated destructor stub
}

void DGSReader::read(std::string path, GraphEventProxy& Gproxy) {

	std::ifstream dgsFile(path);

	std::string line;

	// handle first line
	std::getline(dgsFile, line);
	if (line == "DGS003") {
		DEBUG("found magic cookie: DGS003");
	} else {
		throw std::runtime_error("This does not seem to be a valid DGS file. Expected magic cookie 'DGS003' in first line");
	}

	// handle second line
	std::getline(dgsFile, line);
	std::vector<std::string> split = Aux::StringTools::split(line);

	// TODO: complete implementation

	while (std::getline(dgsFile, line)) {
	}
}

} /* namespace NetworKit */
