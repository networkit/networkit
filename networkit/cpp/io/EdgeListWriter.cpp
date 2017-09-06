/*
 * EdgeListWriter.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "EdgeListWriter.h"
#include "../auxiliary/Log.h"

#include <sstream>

#include "../auxiliary/Enforce.h"

namespace NetworKit {

EdgeListWriter::EdgeListWriter(char separator, node firstNode, bool bothDirections) : separator(separator), firstNode(firstNode), bothDirections(bothDirections) {}

void EdgeListWriter::write(const Graph& G, std::string path) {
	std::ofstream file(path);
	Aux::enforceOpened(file);

	auto writeWeightedEdge = [&](node u, node v, edgeweight weight) {
		file << (u + firstNode) << separator << (v + firstNode) << separator << weight << std::endl;
	};

	auto writeUnweightedEdge = [&](node u, node v){
		file << (u + firstNode) << separator << (v + firstNode) << std::endl;
	};

	if (G.isWeighted()) {
		if (bothDirections) {
			G.forNodes([&](node u) {
				G.forEdgesOf(u, writeWeightedEdge);
			});
		} else {
			G.forEdges(writeWeightedEdge);
		}
	} else {
		if (bothDirections) {
			G.forNodes([&](node u) {
				G.forEdgesOf(u, writeUnweightedEdge);
			});
		} else {
			G.forEdges(writeUnweightedEdge);
		}
	}


	file.close();

}

} /* namespace NetworKit */
