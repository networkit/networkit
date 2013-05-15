/*
 * METISGraphWriter.cpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISGraphWriter.h"

namespace NetworKit {

METISGraphWriter::METISGraphWriter() {
	// TODO Auto-generated constructor stub

}

METISGraphWriter::~METISGraphWriter() {
	// TODO Auto-generated destructor stub
}

void METISGraphWriter::write(Graph& G, std::string path) {
	// TODO: enable weighted graphs

	std::ofstream file(path);
	assert (file.good());

	int64_t n = G.numberOfNodes();
	int64_t m = G.numberOfEdges();

	file << n << " " << m << " " << 0 << std::endl;

	count nc = 0;
	G.forNodes([&](node u) {
		nc += 1;
		G.forNeighborsOf(u, [&](node v){
			file << v << " ";
		});
		file << std::endl;
	});

	assert (nc == n);
	file.close();

}

} /* namespace NetworKit */
