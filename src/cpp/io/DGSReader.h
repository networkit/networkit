/*
 * DGSReader.h
 *
 *  Created on: 01.06.2013
 *      Author: cls
 */

#ifndef DGSREADER_H_
#define DGSREADER_H_

#include <fstream>

#include "DynamicGraphReader.h"

#include <string>
#include <unordered_map>

#include "../dynamics/GraphEventProxy.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

/**
 * @ingroup io
 * DGS is a file format allowing to store graphs and dynamic graphs in a textual human readable way,
 * yet with a small size allowing to store large graphs. Graph dynamics is defined using events like
 * adding, deleting or changing a node or edge. With DGS, graphs will therefore be seen as stream of
 *  such events.
 *
 * Format documentation: http://graphstream-project.org/doc/Advanced-Concepts/The-DGS-File-Format/
 *
 */
class DGSReader: public NetworKit::DynamicGraphReader {

public:
	
	/**
	 * @param[in]	path	Path to file in DGS format.
	 * @param[in]	Gproxy	Graph event proxy receives the events from the file.
	 */
	virtual void read(std::string path, GraphEventProxy& Gproxy);
};

} /* namespace NetworKit */
#endif /* DGSREADER_H_ */
