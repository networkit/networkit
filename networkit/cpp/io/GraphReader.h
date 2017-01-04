/*
 * GraphReader.h
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHREADER_H_
#define GRAPHREADER_H_

#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"


namespace NetworKit {

/**
 * @ingroup io
 * Abstract base class for graph readers.
 */
class GraphReader {
public:
	virtual ~GraphReader() = default;

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	virtual Graph read(const std::string& path) = 0;

        /**
         * Returns whether, based on the given header lines, this reader is able to read the graph file
         * where the header lines where taken from.
         *
         * @param[in] headerLines First lines of a graph file that suffice to uniquely identify the
         *                        graph format being used.
         *                        Example: In the case of an edge list, this should encompass
         *                        1) the line that indicates the number of nodes and edges
         *                        --> Indicates an edge list format in the first place
         *                        2) the first line with with an actual edge definition
         *                        --> Indicates the seperator being used and the lowest node id (usually 0 or 1) 
         *
         * @return whether, based on the given header lines, this reader is able to read the graph file
         */
        virtual bool accepts(const std::vector<std::string&> headerLines) = 0;

};

} /* namespace NetworKit */
#endif /* GRAPHREADER_H_ */
