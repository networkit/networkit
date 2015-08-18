/*
 * DGSStreamParser.cpp
 *
 *  Created on: 23.12.2013
 *      Author: cls
 */

#include "DGSStreamParser.h"
#include "../auxiliary/StringTools.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

DGSStreamParser::DGSStreamParser(std::string path, bool mapped, node baseIndex) : dgsFile(path), mapped(mapped), baseIndex(baseIndex), nextNode(0) {

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
		ERROR("found " , line , " instead of " , cookie , " in first line");
		throw std::runtime_error("expected cookie in first line");
	}

	std::getline(dgsFile, line);
	lc++;
	INFO("DGS stream description: ", line);

	if (mapped) {
		// mapped format

		/**
		 * Maps key string to consecutive, 0-based node id.
		 */
		auto map = [&](std::string key) {
			auto iter = this->key2id.find(key);
			if (iter == key2id.end()) {
				key2id[key] = nextNode;
				nextNode++;
				return key2id[key];
			} else {
				return iter->second;
			}
		};

		while (std::getline(dgsFile, line)) {
			// TRACE(line);
			lc++;
			std::vector<std::string> split = Aux::StringTools::split(line);
			// TRACE("split line: ", split);
			std::string tag = split[0];

			// TODO: remove TRACE

			// parse commands
			if (tag.compare("st") == 0) { // clock
				stream.push_back(GraphEvent(GraphEvent::TIME_STEP));
			} else if (tag.compare("an") == 0) { // add node
				node u = map(split[1]);
				auto ev = GraphEvent(GraphEvent::NODE_ADDITION, u);
				// TRACE(ev.toString());
				stream.push_back(ev);
			} else if (tag.compare("ae") == 0) { // add edge
				node u = map(split[2]);
				node v = map(split[3]);
				edgeweight w = 1.0;
				if (split.size() >= 5) {
					w = std::stod(Aux::StringTools::split(split[4], '=')[1]); // weight=<w>
				}
				auto ev = GraphEvent(GraphEvent::EDGE_ADDITION, u, v, w);
				// TRACE(ev.toString());
				stream.push_back(ev);
			} else if (tag.compare("ce") == 0) { // update edge. Only the "weight" attribute is supported so far
				std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
				node u = map(uvs[0]);
				node v = map(uvs[1]);
				edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
				auto ev = GraphEvent(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, w);
				// TRACE(ev.toString());
				stream.push_back(ev);
			} else if (tag.compare("ie") == 0) { // update edge. Only the "weight" attribute is supported so far
				std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
				node u = map(uvs[0]);
				node v = map(uvs[1]);
				edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
				auto ev = GraphEvent(GraphEvent::EDGE_WEIGHT_INCREMENT, u, v, w);
				// TRACE(ev.toString());
				stream.push_back(ev);
			} else if (tag.compare("de") == 0) {
				std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
				node u = map(uvs[0]);
				node v = map(uvs[1]);
				auto ev = GraphEvent(GraphEvent::EDGE_REMOVAL, u, v);
				// TRACE(ev.toString());
				stream.push_back(ev);
			} else if (tag.compare("dn") == 0) {
				node u = map(split[1]);
				auto ev = GraphEvent(GraphEvent::NODE_REMOVAL, u);
				// TRACE(ev.toString());
				stream.push_back(ev);
			} else if (tag.compare("rn") == 0) {
				node u = map(split[1]);
				auto ev = GraphEvent(GraphEvent::NODE_RESTORATION, u);
				// TRACE(ev.toString());
				stream.push_back(ev);
			}

				else {
				ERROR("malformed line (" , lc , ") : " , line);
				throw std::runtime_error("malformed line in .DGS file");
			}
		}


	} else {
		// direct format
		auto offset = [&](node u) {
			return (u - baseIndex);
		};

		while (std::getline(dgsFile, line)) {
			// TRACE(line);
			lc++;
			std::vector<std::string> split = Aux::StringTools::split(line);
			std::string tag = split[0];

			// parse commands
			if (tag.compare("st") == 0) { // clock
				stream.push_back(GraphEvent(GraphEvent::TIME_STEP));
				// TRACE("read: st ");
			} else if (tag.compare("an") == 0) { // add node
				node u = offset(std::stoul(split[1]));
				stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, u));
				// TRACE("read: an ", u);
			} else if (tag.compare("ae") == 0) { // add edge
				node u = offset(std::stoul(split[2]));
				node v = offset(std::stoul(split[3]));
				edgeweight w = std::stod(Aux::StringTools::split(split[4], '=')[1]); // weight=<w>
				stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, u, v, w));
				// TRACE("read: ae ", u, ",", v, ",", w);
			} else if (tag.compare("ce") == 0) { // update edge. Only the "weight" attribute is supported so far
				std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
				node u = offset(std::stoul(uvs[0]));
				node v = offset(std::stoul(uvs[1]));
				edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
				stream.push_back(GraphEvent(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, w));
				// TRACE("read: ce ", u, ",", v, ",", w);
			} else if (tag.compare("ie") == 0) { // update edge. Only the "weight" attribute is supported so far
				std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
				node u = offset(std::stoul(uvs[0]));
				node v = offset(std::stoul(uvs[1]));
				edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
				stream.push_back(GraphEvent(GraphEvent::EDGE_WEIGHT_INCREMENT, u, v, w));
				// TRACE("read: ce ", u, ",", v, ",", w);
			} else if (tag.compare("de") == 0) {
				std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
				node u = offset(std::stoul(uvs[0]));
				node v = offset(std::stoul(uvs[1]));
				stream.push_back(GraphEvent(GraphEvent::EDGE_REMOVAL, u, v));
				// TRACE("read: de ", u, ",", v);
			} else if (tag.compare("dn") == 0) {
				node u = offset(std::stoul(split[1]));
				stream.push_back(GraphEvent(GraphEvent::NODE_REMOVAL, u));
				// TRACE("read: dn ", u);
			} else if (tag.compare("rn") == 0) {
				node u = offset(std::stoul(split[1]));
				stream.push_back(GraphEvent(GraphEvent::NODE_RESTORATION, u));
				TRACE("read: rn ", u);
			}
			  else {
				ERROR("malformed line (" , lc , ") : " , line);
				throw std::runtime_error("malformed line in .DGS file");
			}
		}

	}




	return stream;
}


} /* namespace NetworKit */
