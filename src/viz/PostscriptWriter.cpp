/*
 * PostscriptWriter.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include "PostscriptWriter.h"

namespace NetworKit {

PostscriptWriter::PostscriptWriter(const Graph& graph, bool isTorus) :
		g(graph), wrapAround(isTorus) {
	numColors = 24;
	psColor = { {1.0, 0.0, 0.0},
		{	1.0, 0.5, 0.0}, {1.0, 1.0, 0.0}, {0.5, 1.0, 0.0}, {0.0, 1.0,
			0.0}, {0.0, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.5, 1.0},
		{	0.0, 0.0, 1.0}, {0.5, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 0.0,
			0.5},

		{	0.6, 0.0, 0.0}, {0.6, 0.3, 0.0}, {0.6, 0.6, 0.0}, {0.3, 0.6,
			0.0}, {0.0, 0.6, 0.0}, {0.0, 0.6, 0.3}, {0.0, 0.6, 0.6},
		{	0.0, 0.3, 0.6}, {0.0, 0.0, 0.6}, {0.3, 0.0, 0.6}, {0.6, 0.0,
			0.6}, {0.6, 0.0, 0.3}};

	ps_minx = g.minCoordinate(0);
	ps_maxx = g.maxCoordinate(0);
	ps_miny = g.minCoordinate(1);
	ps_maxy = g.maxCoordinate(1);

	ps_scale = (ps_sizex - 2 * ps_borderx) / (ps_maxx - ps_minx);
	ps_sizey = (ps_maxy - ps_miny) * ps_scale + 2 * ps_bordery;
}

PostscriptWriter::~PostscriptWriter() {

}

void PostscriptWriter::writeHeader(std::ofstream& file) {
	/* Header */
	if (wrapAround) {
		file << "%!PS-Adobe-3.0 EPSF-3.0\n";
	} else {
		file << "%!PS-Adobe-1.0\n";
	}
	file << "%%Title: clusteredGraph.eps\n";
//	if (wrapAround) {
	file << "%%BoundingBox: 0.000 0.000 1020.0 1020.0\n";
//	}
//	else {
//		file << "%%BoundingBox: 0 0 " << (int)ceil(ps_sizex) << " " << (int)ceil(ps_sizey+ 40.0) << "\n";
//	}
	file << "%%EndComments\n";
	if (!wrapAround) {
		file << "%%EndProlog\n";
		file << "gsave\n";
	}
}

void PostscriptWriter::writeMacros(std::ofstream& file) {
	/* Macros */
	file << "/p {newpath} bind def\n";
	file << "/m {moveto} bind def\n";
	file << "/r {rmoveto} bind def\n";
	file << "/k {rlineto} bind def\n";
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
void PostscriptWriter::writeClustering(Partition& clustering,
		std::ofstream& file) {
	TRACE("start ps writeClustering");

	auto adjust([&](float& val) {
		val += 10.0;
	});

	TRACE("start edge loop in writeClustering");

	/* Kanten zeichnen */
	g.forEdges([&](node u, node v) {
		if (wrapAround) {
			if (clustering[u] == clustering[v] && clustering[u] != none) {
				// same cluster
			float r = psColor[clustering[u] % numColors].r;
			float g = psColor[clustering[u] % numColors].g;
			float b = psColor[clustering[u] % numColors].b;
			file << r << " " << g << " " << b << " c ";
		}
		else {
			// different clusters -> grey
			file << "0.80 0.80 0.80 c 1.0 w ";
		}

		Point<float> start = g.getCoordinate(u);
		Point<float> end = g.getCoordinate(v);
		float& startx = start[0];
		float& starty = start[1];
		float& endx = end[0];
		float& endy = end[1];

		if (wrapAround) {
			auto adjust1([&](float& val) {
						if (val > 500.0f) {
							val -= 1000.0f;
						}
						else if (val < -500.0f) {
							val += 1000.0f;
						}
					});

			auto adjustWrapAround([&](float& distx, float& disty) {
						adjust1(distx);
						adjust1(disty);
					});

			float distx = endx - startx;
			float disty = endy - starty;

			adjustWrapAround(distx, disty);

			endx = startx + distx;
			endy = starty + disty;
		}

		adjust(startx);
		adjust(endx);
		adjust(starty);
		adjust(endy);

		file << "p " << startx << " " << starty << " m " << endx << " " << endy << " l s\n";
//			}
//			else {
//				startx = (startx - ps_minx) * ps_scale + ps_borderx;
//				starty = (starty - ps_miny) * ps_scale + ps_bordery;
//				endx = (endx - ps_minx) * ps_scale + ps_borderx;
//				endy = (endy - ps_miny) * ps_scale + ps_bordery;
//				file << "p " << startx << " " << starty << " m " << endx << " " << endy << " l s\n";
//			}

	}
});

	/* Knoten zeichnen */
	float dotsize = 2.0;
	g.forNodes([&](node u) {
		if (clustering[u] != none) {
			// change color
			float r = psColor[clustering[u] % numColors].r;
			float g = psColor[clustering[u] % numColors].g;
			float b = psColor[clustering[u] % numColors].b;
			file << r << " " << g << " " << b << " c ";
		}
		else {
			file << "0.0 0.0 0.0 c ";
		}

		Point<float> point = g.getCoordinate(u);
		float& x = point[0];
		float& y = point[1];

//		if (wrapAround) {
		adjust(x);
		adjust(y);
		file << "p " << x << " " << y << " " << dotsize << " 0.00 360.00 a s\n";
//		}
//		else {
//			x = (x - ps_minx) * ps_scale + ps_borderx;
// FIXME			y = (y - ps_miny) * ps_scale + ps_bordery;
//			file << "p " << x << " " << y << " " << dotsize << " 0.00 360.00 a\n";
//		}

//		DEBUG("x: " , x);
//		DEBUG("y: " , y);

	});
}

void PostscriptWriter::init(std::string path, std::ofstream& file) {
	TRACE("start ps init");

	file.open(path.c_str());
	if (true || wrapAround) { // FIXME
		file.precision(3);
		// adjust coordinates for postscript output
		g.forNodes([&](node u) {
			TRACE("change coordinate for node " , u);
			g.getCoordinate(u).scale(1000.0);
		});
	} else {
		file.precision(2);
	}
	file << std::fixed;
}

void PostscriptWriter::write(Partition& clustering, std::string path) {
	TRACE("start ps write clustering");

	std::ofstream file;
	init(path, file);
	writeHeader(file);
	writeMacros(file);

	file << "0.00 0.00 0.00 c\n";   // 0.25 0.25 0.25 c 1.0 w

	writeClustering(clustering, file);
	if (!wrapAround) {
		file << "grestore\n";
	}

	file.close();
}

void PostscriptWriter::write(std::string path) {
	TRACE("start ps write");

	std::ofstream file;
	init(path, file);
	writeHeader(file);
	writeMacros(file);

	ClusteringGenerator gen;
	Partition allNone = gen.makeOneClustering(g);
	writeClustering(allNone, file);
	if (!wrapAround) {
		file << "grestore\n";
	}

	file.close();
}


} /* namespace NetworKit */
