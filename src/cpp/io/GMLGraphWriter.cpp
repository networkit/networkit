/*
 * GMLGraphWriter.cpp
 *
 *  Created on: 05.08.2013
 *      Author: Stefan Bertsch
 */

#include "GMLGraphWriter.h"
#include "../auxiliary/Enforce.h"

namespace NetworKit {

void GMLGraphWriter::write(Graph& G, const std::string& path) {
	std::ofstream file(path);
	Aux::enforceOpened(file);
	
	file << "graph\n[\n";
	
	G.forNodes([&](node u) {
		file << "  node\n"
		        "  [\n"
		        "    id " << u << "\n"
		        "  ]\n";
		
		G.forNeighborsOf(u, [&](node v) {
			file << "  edge\n"
			        "  [\n"
			        "    source "<< u << "\n"
			        "    target "<< v << "\n"
			        "  ]\n";
		});
	});
	file << "]\n";
}

} /* namespace NetworKit */
