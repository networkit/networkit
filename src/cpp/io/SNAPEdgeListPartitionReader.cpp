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

Partition SNAPEdgeListPartitionReader::read(std::string path) {
	std::ifstream file(path);

	// check if file readable
	if (!file) {
		throw std::runtime_error("invalid clustering file");
	}

	Partition communities;
	std::string line;
	//	std::stringstream linestream;
	count lines = 1;
	node current;
/*
	while (std::getline(file,line)) {
		std::stringstream linestream(line);
		while (linestream >> current) {
			DEBUG("read nodeID: ", current);
			count currentNElem = communities.numberOfElements();
			DEBUG("currentNElem: ", currentNElem);
			if ( current-1 > currentNElem ) {
				count diff = current - currentNElem-1;
				DEBUG("partition is not big enough; nodeID: ",current, "; currentNElem: ",currentNElem, " and diff: ",diff);
				while (diff >= 0) {
					communities.extend();
					--diff;
				}
				DEBUG("partition size is now: ",communities.numberOfElements());
			}
			communities[current] = lines;
		}
		lines++;
	}
	communities.setUpperBound(lines);*/
	return communities;
}


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


} /* namespace NetworKit */
