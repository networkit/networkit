/*
 * VNAGraphWriter.cpp
 *
 *  Created on: 23.10.2013
 *      Author: Stefan Bertsch
 */

#include "VNAGraphWriter.h"

namespace NetworKit {

VNAGraphWriter::VNAGraphWriter() {

}

VNAGraphWriter::~VNAGraphWriter() {

}

void VNAGraphWriter::write(Graph& G, std::string path) {
	this->write(G, G.isMarkedAsWeighted(), path);
}

void VNAGraphWriter::writeGeneric(Graph& G, bool weighted, std::string path, Partition& partition, count dim) {

	std::ofstream file(path);
	assert (file.good());

//	int64_t n = G.numberOfNodes();
//	int64_t m = G.numberOfEdges();
//	file << n << " " << m << " " << (int)weighted << std::endl;

	file << "*Node properties" << std::endl;
	switch (dim) {
	case 0: { // clustering is defined
		file << "ID x y color" << std::endl;
		break;
	}
	case 2: {
		file << "ID x y" << std::endl;
		break;
	}
	case 3: {
		file << "ID x y z" << std::endl;
		break;
	}
	default: {
		throw std::runtime_error("Dimension not supported by file format GDF. Skip writing file.");
	}
	}

	Point<float> point;
	if (dim == 0) {
		G.forNodes([&](node u) {
			point = G.getCoordinate(u);
			file << u << " " << point[0] << " " << point[1] << " "
				 << 70 * partition.subsetOf(u) << std::endl;
		});
	}
	else {
		G.forNodes([&](node u) {
			point = G.getCoordinate(u);
		    file << u << " " << point.toSsvString();
		});
	}

	file << "*tie data" << std::endl;
	file << "*from to" << std::endl;
	
	G.forEdges([&](node u, node v) {
       	    file << u << " " << v << std::endl;
       	});

	file.close();
}

void VNAGraphWriter::write(Graph& G, bool weighted, std::string path) {
	Partition dummy(0);
	writeGeneric(G, weighted, path, dummy, 2);
}

void VNAGraphWriter::write(Graph& G, bool weighted, std::string path, Partition& partition) {
	writeGeneric(G, weighted, path, partition, 2);
}

void VNAGraphWriter::write3D(Graph& G, bool weighted, std::string path) {
	Partition dummy(0);
	writeGeneric(G, weighted, path, dummy, 3);
}

} /* namespace NetworKit */
