/*
 * FastMETISParser.h
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef FASTMETISPARSER_H_
#define FASTMETISPARSER_H_

#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

class FastMETISParser {
public:
	FastMETISParser();
	virtual ~FastMETISParser();

	Graph parse(const std::string& path);
};

} /* namespace NetworKit */
#endif /* FASTMETISPARSER_H_ */
