/*
 * PostscriptWriter.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include "PostscriptWriter.h"

namespace NetworKit {

PostscriptWriter::PostscriptWriter(Graph& graph): g(&graph) {
	numColors = 24;
	psColor = { { 1.0, 0.0, 0.0 },
			{ 1.0, 0.5, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.5, 1.0, 0.0 }, { 0.0, 1.0,
					0.0 }, { 0.0, 1.0, 0.5 }, { 0.0, 1.0, 1.0 }, { 0.0, 0.5, 1.0 },
			{ 0.0, 0.0, 1.0 }, { 0.5, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 0.0,
					0.5 },

			{ 0.6, 0.0, 0.0 }, { 0.6, 0.3, 0.0 }, { 0.6, 0.6, 0.0 }, { 0.3, 0.6,
					0.0 }, { 0.0, 0.6, 0.0 }, { 0.0, 0.6, 0.3 }, { 0.0, 0.6, 0.6 },
			{ 0.0, 0.3, 0.6 }, { 0.0, 0.0, 0.6 }, { 0.3, 0.0, 0.6 }, { 0.6, 0.0,
					0.6 }, { 0.6, 0.0, 0.3 } };


	ps_minx = g->coordinates.minCoordinate(0);
	ps_maxx = g->coordinates.maxCoordinate(0);
	ps_miny = g->coordinates.minCoordinate(1);
	ps_maxy = g->coordinates.maxCoordinate(1);

	ps_scale = (ps_sizex - 2 * ps_borderx) / (ps_maxx - ps_minx);
	ps_sizey = (ps_maxy - ps_miny) * ps_scale + 2 * ps_bordery;
}

PostscriptWriter::~PostscriptWriter() {
	// TODO Auto-generated destructor stub
}

void PostscriptWriter::writeHeader(std::ofstream& file) {
	/* Header */
	file << "%!PS-Adobe-1.0\n";
	file << "%%Title: clusteredGraph.eps\n";
	file << "%%BoundingBox: 0 0 " << (int)ceil(ps_sizex) << " " << (int)ceil(ps_sizey+ 40.0) << "\n";
	file << "%%EndComments\n";
	file << "%%EndProlog\n";
	file << "gsave\n";
}

void PostscriptWriter::writeMacros(std::ofstream& file) {
	/* Macros */
	file << "/p {newpath} bind def\n";
	file << "/m {moveto} bind def\n";
	file << "/r {rmoveto} bind def\n";
	file << "/l {lineto} bind def\n";
	file << "/n {rlineto} bind def\n";
	file << "/c {setrgbcolor} bind def\n";
	file << "/s {stroke} bind def\n";
	file << "/w {setlinewidth} bind def\n";
	file << "/h {show} bind def\n";
	file << "/a {arc closepath fill} bind def\n";
	file << "/b {closepath eofill} bind def\n";
}

// TODO: node and edge weights and thicker nodes/edges
void PostscriptWriter::writeClustering(Clustering& clustering, std::ofstream& file) {
	/* Kanten zeichnen */
	g->forEdges([&](node u, node v) {
		if (u < v) {
			if (clustering[u] == clustering[v] && clustering[u] != none) {
				// same cluster
				float r = psColor[clustering[u] % numColors].r;
				float g = psColor[clustering[u] % numColors].g;
				float b = psColor[clustering[u] % numColors].b;
				file << r << " " << g << " " << b << " c ";
			}
			else {
				// different clusters -> grey
				file << "0.50 0.50 0.50 c ";
			}

			float startx = (g->coordinates.getCoordinate(u, 0) - ps_minx) * ps_scale + ps_borderx;
			float starty = (g->coordinates.getCoordinate(u, 1) - ps_miny) * ps_scale + ps_bordery;
			float endx = (g->coordinates.getCoordinate(v, 0) - ps_minx) * ps_scale + ps_borderx;
			float endy = (g->coordinates.getCoordinate(v, 1) - ps_miny) * ps_scale + ps_bordery;
			file << "p " << startx << " " << starty << " m " << endx << " " << endy << " l s\n";
		}
	});


	/* Knoten zeichnen */
	float dotsize = 1.2;
	g->forNodes([&](node u) {
		if (clustering[u] != none) {
			// change color
			float r = psColor[clustering[u] % numColors].r;
			float g = psColor[clustering[u] % numColors].g;
			float b = psColor[clustering[u] % numColors].b;
			file << r << " " << g << " " << b << " c ";
		}

		float x = (g->coordinates.getCoordinate(u, 0) - ps_minx) * ps_scale + ps_borderx;
		float y = (g->coordinates.getCoordinate(u, 1) - ps_miny) * ps_scale + ps_bordery;
		file << "p " << x << " " << y << " " << dotsize << " 0.00 360.00 a\n";

	});
}

void PostscriptWriter::write(Clustering& clustering, std::string filename) {
	std::ofstream file;
	file.open(filename.c_str());
	file.precision(2);
	file << std::fixed;

	writeHeader(file);
	writeMacros(file);

	file << "0.00 0.00 0.00 c\n";
	writeClustering(clustering, file);
	file << "grestore\n";

	file.close();
}

} /* namespace NetworKit */
