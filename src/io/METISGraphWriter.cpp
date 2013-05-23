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
	this->write(G, false, path);
}
void METISGraphWriter::write(Graph& G, bool weighted, std::string path) {
	// TODO: enable weighted graphs

	std::ofstream file(path);
	assert (file.good());

	int64_t n = G.numberOfNodes();
	int64_t m = G.numberOfEdges();

	file << n << " " << m << " " << (int)weighted << std::endl;

	if (weighted == false) {
		G.forNodes([&](node u) {
			G.forNeighborsOf(u, [&](node v){
				file << v + 1 << " ";
			});
			file << std::endl;
		});
	} else {
			G.forNodes([&](node u) {
			G.forNeighborsOf(u, [&](node v){
				file << v + 1 << " " << G.weight(u, v) <<"   ";
			});
			file << std::endl;
		});
	}

	file.close();
}

} /* namespace NetworKit */
