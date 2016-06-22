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
#include <map>
#include <fstream>

#include "GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class DGSStreamParser {

public:

	DGSStreamParser(std::string path, bool mapped=true, node baseIndex=0);

	std::vector<GraphEvent> getStream();

private:

	std::ifstream dgsFile;
	bool mapped;
	std::map<std::string, node> key2id;
	node baseIndex;
	node nextNode;


};

} /* namespace NetworKit */

#endif /* DGSSTREAMPARSER_H_ */
