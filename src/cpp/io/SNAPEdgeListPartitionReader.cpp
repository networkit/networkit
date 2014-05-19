/*
 * SNAPEdgeListPartitionReader.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "SNAPEdgeListPartitionReader.h"

namespace NetworKit {

SNAPEdgeListPartitionReader::SNAPEdgeListPartitionReader() {
	// TODO Auto-generated constructor stub

}

SNAPEdgeListPartitionReader::~SNAPEdgeListPartitionReader() {
	// TODO Auto-generated destructor stub
}

Cover SNAPEdgeListPartitionReader::read(std::string path, std::unordered_map<node,node>& mapNodeIds, Graph& G) {
	std::ifstream file;
	std::string line; // the current line

	// read file once to get to the last line and figure out the number of nodes
	// unfortunately there is an empty line at the ending of the file, so we need to get the line before that
	
	file.open(path);
	if (! file.good()) {
		throw std::runtime_error("unable to read from file");
	}
	
//	std::string previousLine;
	node current;
	
	std::string commentPrefix = "#";
	
//	count firstNode = 0;
//	char separator = '\t';

	Cover communities(G.numberOfNodes());

	//DEBUG("separator: " , separator);
	//DEBUG("first node: " , firstNode);

	// first find out the maximum node id
	//DEBUG("first pass");
	count i = 0;
	while (file.good()) {
		++i;
		std::getline(file, line);
		// TRACE("read line: " , line);
		if (line.compare(0, commentPrefix.length(), commentPrefix) == 0) {
			// TRACE("ignoring comment: " , line);
	        } else if (line.length() == 0) {
        		// TRACE("ignoring empty line");
		} else {
			std::stringstream linestream(line);
			while (linestream >> current) {
				if (mapNodeIds.find(current) != mapNodeIds.end()) {
					communities.addToSubset(i,mapNodeIds[current]);
				} else {
					WARN("unknown node ID found (",current,") and ignored");
				}
			}
		}
	}

	file.close();
	communities.setUpperBound(i);
	return communities;
}

#if 0
Partition SNAPEdgeListPartitionReader::readWithInfo(std::string path, count nNodes) {
	std::ifstream file(path);

	// check if file readable
	if (!file) {
		throw std::runtime_error("invalid clustering file");
	}

	Partition communities(nNodes);
	std::string line;
	//	std::stringstream linestream;
	count lines = 1;
	count nValues = 0;
	std::unordered_map<node,count> ids;
	node current;
	while (std::getline(file,line)) {
		std::stringstream linestream(line);
		while (linestream >> current) {
			nValues++;
			if (!ids.insert(std::make_pair(current,0)).second) {
				ids[current]++;
			}
			//DEBUG("read nodeID: ", current);
			if ( current < nNodes ) {
				communities[current] = lines;
			} else {
				DEBUG("found invalid node ID: ",current);
			}
		}
		lines++;
	}
	count sum = 0;
	for(auto x: ids) {
		if (x.second > 1) {
			DEBUG(x.first," has been found ",x.second, " times");
			sum += x.second;
		}
	}
	WARN("found ",ids.size()," values and ",sum," total occurrences");

	communities.setUpperBound(lines);
	return communities;
}
#endif

} /* namespace NetworKit */
