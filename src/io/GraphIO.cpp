/*
 * GraphIO.cpp
 *
 *  Created on: 09.01.2013
 *      Author: cls
 */

#include "GraphIO.h"

namespace EnsembleClustering {

GraphIO::GraphIO() {
	// TODO Auto-generated constructor stub

}

GraphIO::~GraphIO() {
	// TODO Auto-generated destructor stub
}

void GraphIO::toEdgeList(Graph& G, std::string path) {

	std::ofstream file;
	file.open(path.c_str());

	G.forallEdges([&](node v, node w) {
		file << v << " " << w << std::endl;
	});

	file.close();
	INFO("wrote graph to file: " << path);

}

} /* namespace EnsembleClustering */
