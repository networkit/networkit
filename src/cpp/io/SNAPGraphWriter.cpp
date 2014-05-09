/*
 * SNAPGraphWriter.cpp
 *
 *  Created on: 24.09.2013
 *      Author: cls
 */

#include "SNAPGraphWriter.h"

namespace NetworKit {

SNAPGraphWriter::SNAPGraphWriter() {
	// TODO Auto-generated constructor stub

}

SNAPGraphWriter::~SNAPGraphWriter() {
	// TODO Auto-generated destructor stub
}

void SNAPGraphWriter::write(Graph& G, std::string path) {
    std::ofstream file;
    file.open(path);
    assert (file.good());

    // write "problem line" - n, m, directed/undirected, weighted/weight type
    file << "p " << G.numberOfNodes() << " " << G.numberOfEdges() << " u u 0" <<  std::endl; // FIXME: makeshift

    G.forEdges([&](node u, node v){
    	file << u << " " << v << std::endl;
    });

    file.close();
}

} /* namespace NetworKit */
