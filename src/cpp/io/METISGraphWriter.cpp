/*
 * METISGraphWriter.cpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISGraphWriter.h"

#include "../auxiliary/Enforce.h"

namespace NetworKit {

void METISGraphWriter::write(Graph& G, const std::string& path) {
	this->write(G, G.isWeighted(), path);
}
void METISGraphWriter::write(Graph& G, bool weighted, std::string path) {
	std::ofstream file(path);
	Aux::enforceOpened(file);

	int64_t n = G.numberOfNodes();
	int64_t m = G.numberOfEdges();

	file << n << " " << m << " " << int{weighted} << '\n';

	if (weighted == false) {
		G.forNodes([&](node u) {
			G.forNeighborsOf(u, [&](node v){
				file << v + 1 << " ";
			});
			file << '\n';
		});
	} else {
			G.forNodes([&](node u) {
			G.forNeighborsOf(u, [&](node v){
				file << v + 1 << " " << G.weight(u, v) <<"   ";
			});
			file << '\n';
		});
	}
}

} /* namespace NetworKit */
