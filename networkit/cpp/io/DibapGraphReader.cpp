/*
 * DibapGraphReader.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: Henning
 */

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64

#include "DibapGraphReader.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

Graph DibapGraphReader::read(const std::string& path) {
	int n, i;
	short type;
	FILE * file = NULL;

	int V = 0;
	int dvw = 0;
	std::vector<int> vw;
	std::vector<int> of;
	std::vector<int> to;
	int dew = 0;
	std::vector<int> ew;
	int dxy = 0;
	std::vector<float> xy;
	int numE2 = 0;

	// init, try to open file
	DEBUG("reading graph in DibaP format... ");
	file = fopen(path.c_str(), "r");
	if (file == NULL) {
		throw std::runtime_error("cannot open file ");
		return 0;
	}

	// detect file type
	n = fread(&type, sizeof(short), 1, file);
	if (n != 1 || ntohs(type) != IO_TYPE_GI) {
		throw std::runtime_error("bad file structure ");
		return 0;
	}

	// read number of vertices
	n = fread(&V, sizeof(int), 1, file);
	V = ntohl(V);
	TRACE("(V=%d %d) " , V , ", " , n);
	if (n != 1) {
		throw std::runtime_error("bad file structure ");
		return 0;
	}

	// read vertex weight dimension
	n = fread(&dvw, sizeof(int), 1, file);
	dvw = ntohl(dvw);
	TRACE("(dvw=%d %d) " , dvw , ", " , n);

	if (n != 1) {
		throw std::runtime_error("bad file structure ");
		return 0;
	}

	if (dvw > 0) {
		// read vertex weights
		vw.resize(V);
		n = fread(&vw[0], sizeof(int), V * dvw, file);
		TRACE("(vw %d) " , n);
		if (n != V) {
			throw std::runtime_error("bad file structure ");
		}
		for (i = 0; i < n; i++)
			vw[i] = ntohl(vw[i]);
	}

	of.resize(V + 1);
	n = fread(&of[0], sizeof(int), V + 1, file);
	TRACE("(of %d) " , n);
	if (n != V + 1) {
		throw std::runtime_error("bad file structure ");
	}
	for (i = 0; i < n; i++)
		of[i] = ntohl(of[i]);

	numE2 = of[V];
	to.resize(numE2);
	n = fread(&to[0], sizeof(int), numE2, file);
	TRACE("(to %d) " , n);
	if (n != numE2) {
		throw std::runtime_error("bad file structure ");
	}
	for (i = 0; i < n; i++)
		to[i] = ntohl(to[i]);

	n = fread(&dew, sizeof(int), 1, file);
	dew = ntohl(dew);
	TRACE("(dew=%d %d) " , dew , ", " , n);
	if (n != 1) {
		throw std::runtime_error("bad file structure ");
	}

	if (dew > 0) {
		int numWeights = numE2 * dew;
		ew.resize(numWeights);
		n = fread(&ew[0], sizeof(int), numWeights, file);
		TRACE("(ew %d) " , n);
		if (n != numWeights) {
			throw std::runtime_error("bad file structure ");
		}
		for (i = 0; i < n; i++)
			ew[i] = ntohl(ew[i]);
	}

	n = fread(&dxy, sizeof(int), 1, file);
	dxy = ntohl(dxy);
	TRACE("(dxy=%d %d) " , dxy , ", " , n);
	if (n != 1) {
		throw std::runtime_error("bad file structure ");
	}

	if (dxy > 0) {
		int numCoords = V * dxy;
		xy.resize(numCoords);
		n = fread(&xy[0], sizeof(float), numCoords, file);
		TRACE("(xy %d) " , n);
		if (n != numCoords) {
			throw std::runtime_error("bad file structure ");
		}
	}

	fclose(file);

	// fill graph: FIXME: so far without node weights, extend to weights
	Graph graph(0);
	if (dxy == 2) {
		for (index v = 0; v < (count) V; ++v) {
			graph.addNode(xy[2 * v], xy[2 * v + 1]);
		}
	} else {
		for (index v = 0; v < (count) V; ++v) {
			graph.addNode();
		}
	}

	if (dew > 0) {
		for (index v = 0; v < (count) V; ++v) {
			for (index e = of[v]; e < (index) of[v + 1]; ++e) {
				if (v <= (index) to[e]) {
					graph.addEdge(v, to[e], ew[e]);
				}
			}
		}
	} else {
		for (index v = 0; v < (count) V; ++v) {
			for (index e = of[v]; e < (index) of[v + 1]; ++e) {
				if (v <= (index) to[e]) {
					graph.addEdge(v, to[e]);
				}
			}
		}
	}

	graph.shrinkToFit();
	return graph;
}

} /* namespace NetworKit */

#endif // check for non-Windows
