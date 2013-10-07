/*
 * FastMETISParser.h
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef FASTMETISPARSER_H_
#define FASTMETISPARSER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class FastMETISParser {
public:
	FastMETISParser();
	virtual ~FastMETISParser();

	std::vector<std::vector<node>> parse(std::istream& stream);
};

} /* namespace NetworKit */
#endif /* FASTMETISPARSER_H_ */
