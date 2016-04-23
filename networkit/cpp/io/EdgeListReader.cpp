/*
 * EdgeListReader.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "EdgeListReader.h"
#include "../auxiliary/Log.h"

#include <sstream>

#include "../auxiliary/Enforce.h"

namespace NetworKit {

EdgeListReader::EdgeListReader(const char separator, const node firstNode, const std::string commentPrefix, const bool continuous, const bool directed) :
	separator(separator), commentPrefix(commentPrefix), firstNode(firstNode), continuous(continuous), mapNodeIds(), directed(directed) {
//	this->mapNodeIds;i
}

Graph EdgeListReader::read(const std::string& path) {
	if (this->continuous) {
		DEBUG("read graph with continuous ids");
		return readContinuous(path);
	} else {
		DEBUG("read graph with NON continuous ids");
		return readNonContinuous(path);
	}
}

std::map<std::string,node> EdgeListReader::getNodeMap() {
	if (this->continuous) throw std::runtime_error("Input files are assumed to have continuous node ids, therefore no node mapping has been created.");
	return this->mapNodeIds;
}

Graph EdgeListReader::readContinuous(const std::string& path) {
	std::ifstream file(path);
	Aux::enforceOpened(file);
	std::string line; // the current line

	// read file once to get to the last line and figure out the number of nodes
	// unfortunately there is an empty line at the ending of the file, so we need to get the line before that

	node maxNode = 0;
	bool weighted = false;
	bool checkedWeighted = false;

	DEBUG("separator: " , this->separator);
	DEBUG("first node: " , this->firstNode);

	// first find out the maximum node id
	DEBUG("first pass");
	count i = 0;
	while (file.good()) {
		++i;
		std::getline(file, line);
		TRACE("read line: " , line);
		if (!line.empty()) {
			if (line.back() == '\r') line.pop_back();
			if (line.compare(0, this->commentPrefix.length(), this->commentPrefix) == 0) {
				TRACE("ignoring comment: " , line);
			} else {
				std::vector<std::string> split = Aux::StringTools::split(line, this->separator);
				if (!checkedWeighted) {
					if (split.size() == 2) {
						weighted = false;
					} else if (split.size() == 3) {
						INFO("Identified graph as weighted.");
						weighted = true;
					}
					checkedWeighted = true;
				}
				if (split.size() == 2 || split.size() == 3) {
					TRACE("split into : " , split[0] , " and " , split[1]);
					node u = std::stoul(split[0]);
					if (u > maxNode) {
						maxNode = u;
					}
					node v = std::stoul(split[1]);
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
		} else {
			DEBUG("line ", i, " is empty.");
		}
	}
	file.close();
	maxNode = maxNode - this->firstNode + 1;
	DEBUG("max. node id found: " , maxNode);

	Graph G(maxNode, weighted, directed);

	DEBUG("second pass");
	file.open(path);
	// split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes
	i = 0; // count lines
	while(std::getline(file,line)){
		if(*line.rbegin() == '\r') line.pop_back();
		++i;
		if (line.compare(0, this->commentPrefix.length(), this->commentPrefix) == 0) {
			// TRACE("ignoring comment: " , line);
		} else {
			// TRACE("edge line: " , line);
			std::vector<std::string> split = Aux::StringTools::split(line, this->separator);
			if (split.size() == 2) {
				node u = std::stoul(split[0]) - this->firstNode;
				node v = std::stoul(split[1]) - this->firstNode;
			    if (!G.hasEdge(u,v)) {
					G.addEdge(u, v);
				}
			} else if (weighted && split.size() == 3) {
				node u = std::stoul(split[0]) - this->firstNode;
				node v = std::stoul(split[1]) - this->firstNode;
				double weight = std::stod(split[2]);
			    if (!G.hasEdge(u,v)) {
					G.addEdge(u, v, weight);
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

	G.shrinkToFit();
	return G;
}


Graph EdgeListReader::readNonContinuous(const std::string& path) {
	std::ifstream file(path);
	Aux::enforceOpened(file);
	DEBUG("file is opened, proceed");
	std::string line; // the current line
	node consecutiveID = 0;

	bool weighted = false;
	bool checkedWeighted = false;

	// first find out the maximum node id
	DEBUG("first pass: create node ID mapping");
	count i = 0;
	while (file.good()) {
		++i;
		std::getline(file, line);
		TRACE("read line: " , line);
		if (!line.empty()) {
			if(line.back() == '\r') line.pop_back();
			if (line.compare(0, this->commentPrefix.length(), this->commentPrefix) == 0) {
				 TRACE("ignoring comment: " , line);
			} else if (line.length() == 0) {
				TRACE("ignoring empty line");
			} else {
				std::vector<std::string> split = Aux::StringTools::split(line, this->separator);
				if (!checkedWeighted) {
					if (split.size() == 2) {
						weighted = false;
					} else if (split.size() == 3) {
						INFO("Identified graph as weighted.");
						weighted = true;
					}
					checkedWeighted = true;
				}
				if (split.size() == 2 || split.size() == 3) {
					TRACE("split into : " , split[0] , " and " , split[1]);
					if(this->mapNodeIds.insert(std::make_pair(split[0],consecutiveID)).second) ++consecutiveID;
					if(this->mapNodeIds.insert(std::make_pair(split[1],consecutiveID)).second) ++consecutiveID;
				} else {
					std::stringstream message;
					message << "malformed line ";
					message << i << ": ";
					message << line;
					throw std::runtime_error(message.str());
				}
			}
		} else {
			DEBUG("line ", i, " is empty.");
		}
	}
	file.close();

	DEBUG("found ",this->mapNodeIds.size()," unique node ids");
	Graph G(this->mapNodeIds.size(), weighted, directed);

	DEBUG("second pass: add edges");
	file.open(path);

	// split the line into start and end node. since the edges are sorted, the start node has the highest id of all nodes
	i = 0; // count lines
	while(std::getline(file,line)){
		if(*line.rbegin() == '\r') line.pop_back();
        ++i;
		if (line.compare(0, this->commentPrefix.length(), this->commentPrefix) == 0) {
			// TRACE("ignoring comment: " , line);
		} else {
			// TRACE("edge line: " , line);
			std::vector<std::string> split = Aux::StringTools::split(line, this->separator);
			if (split.size() == 2) {
				node u = this->mapNodeIds[split[0]];
				node v = this->mapNodeIds[split[1]];
				if (!G.hasEdge(u,v)) {
					G.addEdge(u, v);
				}
			} else if (weighted && split.size() == 3) {
				node u = this->mapNodeIds[split[0]];
				node v = this->mapNodeIds[split[1]];
				double weight = std::stod(split[2]);
			    if (!G.hasEdge(u,v)) {
					G.addEdge(u, v, weight);
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

} /* namespace NetworKit */
