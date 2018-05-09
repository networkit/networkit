/*
 * FastSNAPGraphReader.cpp
 *
 *  Created on: 04.05.2018
 *      Author: Alexander van der Grinten
 */

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "FastSNAPGraphReader.h"
#include "../auxiliary/StringTools.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

Graph FastSNAPGraphReader::read(const std::string &path) {
	Graph G;

	// Maps SNAP node IDs to consecutive NetworKit node IDs.
	auto mapNode = [&] (node in) -> node {
		auto it = nodeIdMap.find(in);
		if(it != nodeIdMap.end())
			return it->second;
		auto result = nodeIdMap.insert({in, G.addNode()});
		if(!result.second)
			throw std::runtime_error("Error in mapping nodes");
		return result.first->second;
	};

	// This function modifies the graph on input.
	auto handleEdge = [&] (node source, node target) {
		if(!G.hasEdge(source, target))
			G.addEdge(source, target);
	};

	auto fd = open(path.c_str(), O_RDONLY);
	if(fd < 0)
		throw std::runtime_error("Unable to open file");

	struct stat st;
	if(fstat(fd, &st))
		throw std::runtime_error("Could not obtain file stats");

	// It does not really matter if we use a private or shared mapping.
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

		//if(!(*it == '\n')){
		//		std::string str(it);
		//		throw std::runtime_error("Error in parsing file");
		//}
		++it;
	}

	if(munmap(window, st.st_size))
		throw std::runtime_error("Could not unmap file");

	G.shrinkToFit();
	return G;
}

} // namespace NetworKit
