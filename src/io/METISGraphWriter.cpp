/*
 * METISGraphWriter.cpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISGraphWriter.h"

namespace EnsembleClustering {

METISGraphWriter::METISGraphWriter() {
	// TODO Auto-generated constructor stub

}

METISGraphWriter::~METISGraphWriter() {
	// TODO Auto-generated destructor stub
}

void METISGraphWriter::write(Graph& G, std::string path) {

	std::ofstream file(path);
	assert (file.good());

	int64_t n = G.numberOfNodes();
	int64_t m = G.numberOfEdges();

	file << n << " " << m << std::endl;

	G.forNodes([&](node u) {
		G.forNeighborsOf(u, [&](node v){
			file << v << " ";
		});
		file << std::endl;
	});

	file.close();

}

} /* namespace EnsembleClustering */
