/*
 * DGSStreamParser.h
 *
 *  Created on: 23.12.2013
 *      Author: cls
 */

#ifndef DGSSTREAMPARSER_H_
#define DGSSTREAMPARSER_H_

#include <string>
#include <vector>
#include <fstream>

#include "GraphEvent.h"


namespace NetworKit {

class DGSStreamParser {

public:

	DGSStreamParser(std::string path);

	virtual std::vector<GraphEvent> getStream();

private:

	std::ifstream dgsFile;

};

} /* namespace NetworKit */

#endif /* DGSSTREAMPARSER_H_ */
