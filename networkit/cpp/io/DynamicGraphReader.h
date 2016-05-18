/*
 * DynamicGraphReader.h
 *
 *  Created on: 01.06.2013
 *      Author: cls
 */

#ifndef DYNAMICGRAPHREADER_H_
#define DYNAMICGRAPHREADER_H_

#include "../dynamics/GraphEventProxy.h"

namespace NetworKit {

/**
 * @ingroup io
 */
class DynamicGraphReader {

public:
	/**
	 * @param[in]	path	path to dynamic graph file
	 * @param[in]	Gproxy	graph event proxy receives the events from the file
	 */
	virtual void read(std::string path, GraphEventProxy& Gproxy) =0;
};

} /* namespace NetworKit */
#endif /* DYNAMICGRAPHREADER_H_ */
