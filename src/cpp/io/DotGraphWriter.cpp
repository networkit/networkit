/*
 * DotWriter.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#include "DotGraphWriter.h"

namespace NetworKit {

DotGraphWriter::DotGraphWriter() {
	// TODO Auto-generated constructor stub

}

DotGraphWriter::~DotGraphWriter() {
	// TODO Auto-generated destructor stub
}

void DotGraphWriter::write(Graph& graph, std::string path) const {
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
