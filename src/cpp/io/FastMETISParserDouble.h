/*
 * FastMETISParser.h
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef FASTMETISPARSERDOUBLE_H_
#define FASTMETISPARSERDOUBLE_H_

#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

class FastMETISParserDouble {
public:
	FastMETISParserDouble();
	virtual ~FastMETISParserDouble();

	Graph parse(const std::string& path);
};

} /* namespace NetworKit */
#endif /* FASTMETISPARSER_H_ */
