/*
 * GMLGraphWriter.cpp
 *
 *  Created on: 05.08.2013
 *      Author: Stefan Bertsch
 */

#include "GMLGraphWriter.h"

namespace NetworKit {

GMLGraphWriter::GMLGraphWriter() {

}

GMLGraphWriter::~GMLGraphWriter() {

}

void GMLGraphWriter::write(Graph& G, std::string path) {
	std::ofstream file(path);
	assert (file.good());

	file << "graph" << std::endl;
	file << "["<< std::endl;

	G.forNodes([&](node u) {
	    file << "  node" << std::endl;
	    file << "  [" << std::endl;
	    file << "    id "<< u << std::endl;
		file << "  ]"<< std::endl;

		G.forNeighborsOf(u, [&](node v) {
		    file <<"  edge" << std::endl;
       	    file << "  [" << std::endl;
       	    file << "    source "<< u << std::endl;
			file << "    target "<< v << std::endl;
			file << "  ]" << std::endl;
		});
	});

	file << "]" << std::endl;
	file.close();
}

} /* namespace NetworKit */
