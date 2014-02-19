/*
 * ClusteringWriter.cpp
 *
 *  Created on: 22.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringWriter.h"

namespace NetworKit {

ClusteringWriter::ClusteringWriter() {
	// TODO Auto-generated constructor stub

}

ClusteringWriter::~ClusteringWriter() {
	// TODO Auto-generated destructor stub
}

void ClusteringWriter::write(Partition& zeta, std::string path) const {
	std::ofstream file;
	file.open(path.c_str());

	zeta.forEntries([&](node v, index c){
		file << c << std::endl;
	});

	file.close();
}

} /* namespace NetworKit */
