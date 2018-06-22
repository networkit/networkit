/*
 * KONECTGraphReader.cpp
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 *
 * The KONECT format is described in detail in
 * http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
 *
 */

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include "../auxiliary/Log.h"

#include "KONECTGraphReader.h"

namespace NetworKit{
	KONECTGraphReader::KONECTGraphReader(bool remapNodes, MultipleEdgesHandling handlingmethod):
	remapNodes(remapNodes), multipleEdgesHandlingMethod(handlingmethod){}

	Graph KONECTGraphReader::read(const std::string &path){
		std::string graphFormat = "";
		std::string graphType = "";
		count numberOfNodes = -1;
		count numberOfEdges = -1;

		bool directed = true;
		bool weighted = false;
		bool multiple = false;
		int valuesPerLine = 2;
		bool secondPropertyLine = false;
		std::unordered_map<node, node> nodeIdMap;

		//open file
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
			if(it == end){ // proper error message if file ends unexpected
				throw std::runtime_error("Unexpected end of file. Whitespace at end of file found");
			}
		};

		// This function parses in whole words
		auto scanWord = [&] () {
			std::string word = "";
			while(it != end && *it != ' ' && *it != '\t' && *it != '\n'){
				word += *it;
				++it;
			}
			return word;
		};

		auto scanId = [&] () -> node {
			char *past;
			auto value = strtol(it, &past, 10);
			if(past <= it){
				throw std::runtime_error("Scanning node failed. The file may be corrupt.");
			}
			it = past;
			return value;
		};

		auto scanWeight = [&] () -> edgeweight{
			char *past;
			auto value = strtod(it, &past);
			if(past <= it){
				throw std::runtime_error("Scanning node failed. The file may be corrupt.");
			}
			it = past;
			return value;
		};

		//parse graph properties
		if(*it == '%'){
			++it;
			skipWhitespace();

			// Parse graph format directed / undirected
			graphFormat = scanWord();

			if (graphFormat == "sym" || graphFormat == "bip"){
				directed = false;
			}else if(graphFormat != "asym"){
				throw std::runtime_error("Graph format not tagged properly. Format is: " + graphFormat);
			}
			skipWhitespace();
			// Parse graph weighting
			graphType = scanWord();
			//NOTE: This parser disregards the edge weight classification of "interval scale" and "ratio scale"
			if (graphType == "weighted" || graphType == "posweighted" || graphType == "signed"){
				weighted = true;
				valuesPerLine = 3;
			} else if (graphType == "positive"){ //multiple edges
				if(multipleEdgesHandlingMethod == SUM_WEIGHTS_UP){
					weighted = true;
				}
				multiple = true; //weights will be added
			} else if (graphType == "multisigned" || graphType == "multiweighted" || graphType == "multiposweighted"){
				weighted = true;
				multiple = true;
				valuesPerLine = 3;
			} else if (graphType == "dynamic"){
				throw std::runtime_error("Dynamic graphs are not supported yet");
			} else if (graphType != "unweighted"){
				throw std::runtime_error("Graph weight not tagged properly. Weight is: " + graphType);
			}
			DEBUG("First property line read in. Format: "+graphFormat+" / Type: "+graphType);
			if (graphFormat == "bip"){
				INFO("Your graph is flagged as a bipartite one. Keep in mind that"
				 "NetworKit cannot handle this kind of graph in its special cases."
				 "It is imported as a usual undirected graph and edges might be discarded.");
			}
			if(multiple){
				DEBUG("Selected handling method for multiple edges is: "+std::to_string(multipleEdgesHandlingMethod));
			}
		} else {
			throw std::runtime_error("No graph properties line found. This line is mandatory.");
		}

		skipWhitespace();
		if (*it != '\n'){
			throw std::runtime_error("No break symbol after first property line");
		}else{
			++it;
		}
		skipWhitespace();
		//second optional property line
		if(*it == '%'){
			++it;
			skipWhitespace();
			numberOfEdges = scanId();
			skipWhitespace();
			numberOfNodes = scanId();
			secondPropertyLine = true;
			while(it != end && *it != '\n')
				++it;
			if(it >= end){ // proper error message if file ends unexpected
				throw std::runtime_error("Unexpected end of file");
			}
			if (*it != '\n'){
				throw std::runtime_error("File not properly formatted. Last character in second property line is: "+std::string(1, *it));
			}else{
				++it;
			}
			DEBUG("Second property line read in. Edges: "+std::to_string(numberOfEdges)+ " / Nodes: "+std::to_string(numberOfNodes));
		}


		Graph graph((secondPropertyLine ? numberOfNodes : 0), weighted, directed);

		//Map nodes and increase graph size if no second property is defined
		std::function<node(node)> mapNode = [&] (node in) -> node {
			if(secondPropertyLine){
				if (in > numberOfNodes){ //if file is corrupted
					//secondPropertyLine is made useless
					ERROR("Given number of nodes by file does not match actual graph");
					secondPropertyLine = false;
					if(remapNodes){ // if remapNodes is selected true, the map has to be initalized with the existing nodes
						nodeIdMap.reserve(numberOfNodes);
						graph.forNodes([&](node u) {
							nodeIdMap.insert({u,u});
						});
					}
					return mapNode(in);
				}else{
					return in - 1; //minus firstNode
				}
			}else{
				if(remapNodes){
					auto it = nodeIdMap.find(in);
					if(it != nodeIdMap.end())
						return it->second;
					auto result = nodeIdMap.insert({in, graph.addNode()});
					assert(result.second);
					return result.first->second;
				}else{
					for(count i = graph.upperNodeIdBound(); i < in; i++)
						graph.addNode();
					return in - 1;
				}
			}
		};

		//Helper function for handling edges
		auto handleEdge = [&] (node source, node target, edgeweight weight = defaultEdgeWeight) {
			if(!graph.hasEdge(source,target)){
				graph.addEdge(source, target, weight);
			}else if (multiple){
				switch(multipleEdgesHandlingMethod){
					case DISCARD_EDGES: break; //Do nothing
					case SUM_WEIGHTS_UP:
						graph.increaseWeight(source, target, weight);
						break;
					case KEEP_MINIUM_WEIGHT:
						if(graph.weight(source,target) > weight){
							graph.setWeight(source,target, weight);
						}
						break;
					default:
						throw std::runtime_error("Invalid multipleEdgesHandlingMethod: "
						+std::to_string(multipleEdgesHandlingMethod));
				}
			}else{
				DEBUG("["+std::to_string(source)+"->"+std::to_string(target)+
				"] Multiple edges detected but declared as: "+graphFormat+
				" and "+graphType);
			}
		};

		while(it != end){
			skipWhitespace();
			if(*it == '\n') {
				// We ignore empty lines.
			} else if(*it == '#' || *it == '%') {
				// Skip non-linebreak characters.
				while(it != end && *it != '\n')
					++it;
			} else {
				// Normal case parsing
				auto sourceId = scanId();
				if(it >= end){ // proper error message if file ends unexpected
					throw std::runtime_error("Unexpected end of file");
				}
				if(!(*it == ' ' || *it == '\t')){
					throw std::runtime_error("Target ID cannot be parsed, maybe the file is corrupt");
				}
				skipWhitespace();

				auto targetId = scanId();
				if(it >= end){ // proper error message if file ends unexpected
					throw std::runtime_error("Unexpected end of file");
				}
				skipWhitespace();
				if (valuesPerLine > 2){
					auto edgeWeight = scanWeight();
					handleEdge(mapNode(sourceId), mapNode(targetId), edgeWeight);
				}else{
					handleEdge(mapNode(sourceId), mapNode(targetId));
				}
			}
			//break lines
			skipWhitespace();
			if (*it != '\n'){
				throw std::runtime_error("File not properly formatted. Last character in line is: "+std::string(1, *it));
			}
			++it;
		}
		//unmap file
		if(munmap(window, st.st_size))
			throw std::runtime_error("Could not unmap file");

		graph.shrinkToFit();
		return graph;
	}
} //namespace NetworKit
