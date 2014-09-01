/*
 * VNAGraphWriter.cpp
 *
 *  Created on: 23.10.2013
 *      Author: Stefan Bertsch
 */

#include "VNAGraphWriter.h"
#include "../auxiliary/Enforce.h"

namespace NetworKit {

void VNAGraphWriter::write(Graph& G, const std::string& path) {
	this->write(G, G.isWeighted(), path);
}

void VNAGraphWriter::writeGeneric(Graph& G, bool weighted, const std::string& path,
                                  Partition& partition, count dim) {

	std::ofstream file(path);
	Aux::enforceOpened(file);

//	int64_t n = G.numberOfNodes();
//	int64_t m = G.numberOfEdges();
//	file << n << " " << m << " " << (int)weighted << '\n';

	file << "*Node properties\n";
	switch (dim) {
	case 0: { // clustering is defined
		file << "ID x y color\n";
		break;
	}
	case 2: {
		file << "ID x y\n";
		break;
	}
	case 3: {
		file << "ID x y z\n";
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
				 << 70 * partition.subsetOf(u) << '\n';
		});
	}
	else {
		G.forNodes([&](node u) {
			point = G.getCoordinate(u);
			file << u << " " << point.toSsvString();
		});
	}

	file << "*tie data\n";
	file << "*from to\n";
	
	G.forEdges([&](node u, node v) {
		file << u << " " << v << '\n';
	});

}

void VNAGraphWriter::write(Graph& G, bool weighted, const std::string& path) {
	Partition dummy(0);
	writeGeneric(G, weighted, path, dummy, 2);
}

void VNAGraphWriter::write(Graph& G, bool weighted, const std::string& path,
                           Partition& partition) {
	writeGeneric(G, weighted, path, partition, 2);
}

void VNAGraphWriter::write3D(Graph& G, bool weighted, const std::string& path) {
	Partition dummy(0);
	writeGeneric(G, weighted, path, dummy, 3);
}

} /* namespace NetworKit */
