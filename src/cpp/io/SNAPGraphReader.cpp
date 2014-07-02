/*
 * SNAPGraphReader.cpp
 *
 *  Created on: 19.05.2014
 *      Author: Maximilian Vogel
 */

#include "SNAPGraphReader.h"
#include "../auxiliary/StringTools.h"
#include "../auxiliary/Log.h"

#include <sstream>
#include <fstream>

namespace NetworKit {

Graph SNAPGraphReader::read(const std::string& path) {
	std::ifstream file;
	std::string line; // the current line

	// read file once to get to the last line and figure out the number of nodes
	// unfortunately there is an empty line at the ending of the file, so we need to get the line before that
	
	file.open(path);
	if (! file.good()) {
		throw std::runtime_error("unable to read from file");
	}
	
	std::string previousLine;
	node maxNode = 0;
	node consecutiveID = 0;
	//std::unordered_map<node,node> mapNodeIds;
	
	std::string commentPrefix = "#";
	
	// count firstNode = 0;
	char separator = '\t';

	//DEBUG("separator: " , separator);
	//DEBUG("first node: " , firstNode);

	// first find out the maximum node id
	DEBUG("first pass: create node ID mapping");
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
			std::vector<std::string> split = Aux::StringTools::split(line, separator);
	
			if (split.size() == 2) {
        			TRACE("split into : " , split[0] , " and " , split[1]);
				node u = std::stoul(split[0]);
				if(mapNodeIds.insert(std::make_pair(u,consecutiveID)).second) consecutiveID++;
				if (u > maxNode) {
					maxNode = u;
				}
				node v = std::stoul(split[1]);
				if(mapNodeIds.insert(std::make_pair(v,consecutiveID)).second) consecutiveID++;

				if (v > maxNode) {
					maxNode = v;
				}
			} else {
				std::stringstream message;
				message << "malformed line ";
				message << i << ": ";
				message << line;
				throw std::runtime_error(message.str());
			}
		}
	}

	file.close();

	//maxNode = maxNode - firstNode + 1;
	//DEBUG("max. node id found: " , maxNode);

	//Graph G(maxNode);
	DEBUG("found ",mapNodeIds.size()," unique node ids");
	Graph G(mapNodeIds.size());

	DEBUG("second pass: add edges");
	file.open(path);

    // split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes
	i = 0; // count lines
	while(std::getline(file,line)){
        	++i;
		if (line.compare(0, commentPrefix.length(), commentPrefix) == 0) {
			// TRACE("ignoring comment: " , line);
		} else {
			// TRACE("edge line: " , line);
			std::vector<std::string> split = Aux::StringTools::split(line, separator);
			std::string splitZero = split[0];
			if (split.size() == 2) {
				node u = mapNodeIds[std::stoi(split[0])];
				node v = mapNodeIds[std::stoi(split[1])];
				if (!G.hasEdge(u,v) && !G.hasEdge(v,u)) {
					G.addEdge(u, v);
				}
			} else {
				std::stringstream message;
				message << "malformed line ";
				message << i << ": ";
				message << line;
				throw std::runtime_error(message.str());
			}
		}
	}
	DEBUG("read ",i," lines and added ",G.numberOfEdges()," edges");
	file.close();

	G.shrinkToFit();
	return G;
}

std::unordered_map<node,node> SNAPGraphReader::getNodeIdMap() {
	return mapNodeIds;
}



}

