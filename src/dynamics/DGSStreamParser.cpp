/*
 * DGSStreamParser.cpp
 *
 *  Created on: 23.12.2013
 *      Author: cls
 */

#include "DGSStreamParser.h"
#include "../auxiliary/StringTools.h"


namespace NetworKit {

DGSStreamParser::DGSStreamParser(std::string path) : dgsFile(path) {

}

std::vector<GraphEvent> NetworKit::DGSStreamParser::getStream() {
	if (! dgsFile.is_open()) {
		throw std::runtime_error("DGS input file could not be opened.");
	}

	std::vector<GraphEvent> stream; // stream containing the events

	std::string line;
	count lc = 0; // line count

	std::string cookie = "DGS004";
	std::getline(dgsFile, line); // get DGS version
	lc++;
	if (line.compare(0, cookie.size(), cookie)) {
		ERROR("found " << line << " instead of " << cookie << " in first line");
		throw std::runtime_error("expected cookie in first line");
	}

	std::getline(dgsFile, line);
	lc++;
	// TODO: handle second line


	while (std::getline(dgsFile, line)) {
		TRACE(line);
		lc++;
		std::vector<std::string> split = Aux::StringTools::split(line);
		std::string tag = split[0];

		// parse commands
		if (tag.compare("st") == 0) { // clock
			stream.push_back(GraphEvent(GraphEvent::TIME_STEP));
		} else if (tag.compare("an") == 0) { // add node
			node u = std::stoul(split[1]);
			stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, u));
		} else if (tag.compare("ae") == 0) { // add edge
			node u = std::stoul(split[2]);
			node v = std::stoul(split[3]);
			edgeweight w = std::stod(Aux::StringTools::split(split[4], '=')[1]); // weight=<w>
			stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, u, v, w));
		} else if (tag.compare("ce") == 0) { // update edge. Only the "weight" attribute is supported so far
			std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
			TRACE(uvs[0] << " " << uvs[1]);
			node u = std::stoul(uvs[0]);
			node v = std::stoul(uvs[1]);
			edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
			stream.push_back(GraphEvent(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, w));
		} else if (tag.compare("de") == 0) {
			std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
			node u = std::stoul(uvs[0]);
			node v = std::stoul(uvs[1]);
			stream.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, u, v));
		} else if (tag.compare("dn") == 0) {
			node u = std::stoul(split[1]);
			stream.push_back(GraphEvent(GraphEvent::NODE_REMOVAL, u));
		} else {
			ERROR("malformed line (" << lc << ") : " << line);
			throw std::runtime_error("malformed line in .DGS file");
		}

	}
	return stream;
}


} /* namespace NetworKit */

