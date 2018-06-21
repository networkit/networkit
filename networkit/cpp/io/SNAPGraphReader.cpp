/*
 * SNAPGraphReader.cpp
 *
 *  Created on: 04.05.2018
 *      Author: Alexander van der Grinten
 */

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "SNAPGraphReader.h"
#include "../auxiliary/StringTools.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

SNAPGraphReader::SNAPGraphReader(bool directed, const bool& remapNodes, const count& nodeCount) :
directed(directed), nodeCount(nodeCount), remapNodes(remapNodes){}

Graph SNAPGraphReader::read(const std::string &path) {
	Graph graph(0,false,directed);

	//In the actual state this parameter has very little influence on the reader performance.
	//There can be a significant boost if it is possible to reserve space in the graph initialization
	if(nodeCount != 0 && remapNodes)
		nodeIdMap.reserve(nodeCount);

	// Maps SNAP node IDs to consecutive NetworKit node IDs.
	auto mapNode = [&] (node in) -> node {
		if (remapNodes){
			auto it = nodeIdMap.find(in);
			if(it != nodeIdMap.end())
				return it->second;
			auto result = nodeIdMap.insert({in, graph.addNode()});
			if(!result.second)
				throw std::runtime_error("Error in mapping nodes");
			return result.first->second;
		}
		for(count i = graph.upperNodeIdBound(); i < in + 1; i++)
			graph.addNode();
		return in;
	};

	// This function modifies the graph on input.
	auto handleEdge = [&] (node source, node target) {
		if(!graph.hasEdge(source, target)){
			graph.addEdge(source, target);
		}else{
			DEBUG("["+std::to_string(source)+"->"+std::to_string(target)+
				"] Multiple edges detected");
		}
	};

	auto fd = open(path.c_str(), O_RDONLY);
	if(fd < 0)
		throw std::runtime_error("Unable to open file");

	struct stat st;
	if(fstat(fd, &st))
		throw std::runtime_error("Could not obtain file stats");

	// It does not really matter if we use a private or shared mappingraph.
	auto window = mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if(window == reinterpret_cast<void *>(-1))
		throw std::runtime_error("Could not map file");

	if(close(fd))
		throw std::runtime_error("Error during close()");

	auto it = reinterpret_cast<char *>(window);
	auto end = reinterpret_cast<char *>(window) + st.st_size;

	// The following functions are helpers for parsing.
	auto skipWhitespace = [&] {
		while(it != end && (*it == ' ' || *it == '\t'))
			++it;
	};

	auto scanId = [&] () -> node {
		char *past;
		auto value = strtol(it, &past, 10);
		if(past <= it)
			throw std::runtime_error("Error in parsing file - looking for nodeId failed");
		it = past;
		return value;
	};

	// This loop does the actual parsing.
	while(it != end) {
		if(it >= end)
			throw std::runtime_error("Unexpected end of file");

		skipWhitespace();

		if(it == end)
			throw std::runtime_error("Unexpected end of file");

		if(*it == '\n') {
			// We ignore empty lines.
		}else if(*it == '#') {
			// Skip non-linebreak characters.
			while(it != end && *it != '\n')
				++it;
		}else{
			auto sourceId = scanId();
			if(it == end)
				throw std::runtime_error("Unexpected end of file");
			if(!(*it == ' ' || *it == '\t'))
				throw std::runtime_error("Error in parsing file - pointer is whitespace");

			skipWhitespace();

			auto targetId = scanId();
			skipWhitespace();

			handleEdge(mapNode(sourceId), mapNode(targetId));
		}

		if(it == end)
			throw std::runtime_error("Unexpected end of file");
		if(!(*it == '\n' || *it == '\r'))
			throw std::runtime_error("Line does not end with line break");
		++it;
	}

	if(munmap(window, st.st_size))
		throw std::runtime_error("Could not unmap file");
	graph.shrinkToFit();
	return graph;
}

} // namespace NetworKit
