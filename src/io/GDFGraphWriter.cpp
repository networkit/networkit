/*
 * GDFGraphWriter.cpp
 *
 *  Created on: 27.10.2013
 *      Author: Stefan Bertsch
 */

#include "GDFGraphWriter.h"

namespace NetworKit {

GDFGraphWriter::GDFGraphWriter() {

}

GDFGraphWriter::~GDFGraphWriter() {

}

void GDFGraphWriter::write(Graph& G, std::string path) {
	this->write(G, G.isMarkedAsWeighted(), path);
}

void GDFGraphWriter::writeGeneric(Graph& G, bool weighted, std::string path, count dim) {
	std::ofstream file(path);
	assert(file.good());

//	int64_t n = G.numberOfNodes();
//	int64_t m = G.numberOfEdges();
//	file << n << " " << m << " " << (int)weighted << std::endl;

	switch (dim) {
	case 2: {
		file << "*nodedef>name VARCHAR,x DOUBLE,y DOUBLE" << std::endl;
		break;
	}
	case 3: {
		file << "nodedef>name VARCHAR,x DOUBLE,y DOUBLE,z DOUBLE" << std::endl;
		break;
	}
	default: {
		WARN("Dimension ", dim, " not supported by file format GDF. Skip writing file.");
		file.close();
		return;
	}
	}

	Point<float> point; // HM: needs to be float since Graph uses float
	G.forNodes([&](node u) {
		point = G.getCoordinate(u);
		file << u << "," << point.toCsvString() << std::endl;
	});

	file << "edgedef>node1 VARCHAR,node2 VARCHAR" << std::endl;

	G.forEdges([&](node u, node v) {
		file << u << "," << v << std::endl;
	});

	file.close();
}

void GDFGraphWriter::write(Graph& G, bool weighted, std::string path) {
	writeGeneric(G, weighted, path, 2);
}

void GDFGraphWriter::write3D(Graph& G, bool weighted, std::string path) {
	writeGeneric(G, weighted, path, 3);
}

} /* namespace NetworKit */
