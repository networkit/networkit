/*
 * DotWriter.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#include "DotWriter.h"

namespace NetworKit {

DotWriter::DotWriter() {
	// TODO Auto-generated constructor stub

}

DotWriter::~DotWriter() {
	// TODO Auto-generated destructor stub
}

void DotWriter::write(Graph& graph, std::string path) const {
	std::ofstream file;
		file.open(path.c_str());

	    file << "graph {" << std::endl;


	    graph.forEdges([&](node u, node v){
	        file << u << " -- " << v << ";" << std::endl;
	    });

	    file << "}" << std::endl;

		file.close();
}

} /* namespace NetworKit */
