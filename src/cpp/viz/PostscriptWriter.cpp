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

	TRACE("ps_minx: ", ps_minx);
	TRACE("ps_maxx: ", ps_maxx);
	TRACE("ps_miny: ", ps_miny);
	TRACE("ps_maxy: ", ps_maxy);
	TRACE("ps_scale: ", ps_scale);
	TRACE("ps_sizey: ", ps_sizey);
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
	file << "%%Title: NetworKit visualization\n";
//	if (wrapAround) {
	file << "%%BoundingBox: 0.000 0.000 1020.0 1020.0\n";
//	}
//	else {
//		file << "%%BoundingBox: 0 0 " << (int)ceil(ps_sizex) << " " << (int)ceil(ps_sizey+ 40.0) << "\n";
//	}
	file << "%%EndComments\n";
	if (! wrapAround) {
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
// FIXME: remove		val += 10.0;
	});

	TRACE("start edge loop in writeClustering, wrapAround? ", wrapAround);
	TRACE("num edges: ", g.numberOfEdges());

	/* Kanten zeichnen */
	g.forEdges([&](node u, node v) {

		if (clustering[u] == clustering[v] && clustering[u] != none) {
			// same cluster
			float r = psColor[clustering[u] % numColors].r;
			float g = psColor[clustering[u] % numColors].g;
			float b = psColor[clustering[u] % numColors].b;
			file << r << " " << g << " " << b << " c ";
//			TRACE("set color to ", r, " ", g, " ", b);
		}
		else {
			// different clusters -> grey
			file << "0.80 0.80 0.80 c 1.0 w ";
//			TRACE("set color to grey");
		}

		Point<float> start = g.getCoordinate(u);
		Point<float> end = g.getCoordinate(v);
		float& startx = start[0];
		float& starty = start[1];
		float& endx = end[0];
		float& endy = end[1];

		auto adjust1([&](float& val) { // TODO: externalize constants
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

		if (wrapAround) {
			adjustWrapAround(distx, disty);
		}

		endx = startx + distx;
		endy = starty + disty;

		adjust(startx);
		adjust(endx);
		adjust(starty);
		adjust(endy);

//		else {
//			startx = (startx - ps_minx) * ps_scale + ps_borderx;
//			starty = (starty - ps_miny) * ps_scale + ps_bordery;
//			endx = (endx - ps_minx) * ps_scale + ps_borderx;
//			endy = (endy - ps_miny) * ps_scale + ps_bordery;
//		}

		file << "p " << startx << " " << starty << " m " << endx << " " << endy << " l s\n";
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
//		TRACE("write coordinate to file: ", x, ", ", y);

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

	float ps_stretchx = (ps_sizex - 2 * ps_borderx);
	float ps_stretchy = (ps_sizey - 2 * ps_bordery);

	auto adjustCoordinate([&](Point<float>& p) {
		p[0] -= ps_minx;
		p[1] -= ps_miny;

		p[0] *=  ps_stretchx / (ps_maxx - ps_minx);
		p[1] *=  ps_stretchy / (ps_maxy - ps_miny);

		p[0] += ps_borderx;
		p[1] += ps_bordery;

		TRACE("New coordinate: ", p.toCsvString());

		return p;
	});

	file.open(path.c_str());
	file.precision(3);
	// adjust coordinates for postscript output
	TRACE("change coordinates for nodes to fit into bounding box");
	g.forNodes([&](node u) {
		g.setCoordinate(u, adjustCoordinate(g.getCoordinate(u)));
	});

	ps_minx = g.minCoordinate(0);
	ps_maxx = g.maxCoordinate(0);
	ps_miny = g.minCoordinate(1);
	ps_maxy = g.maxCoordinate(1);

	file << std::fixed;

	writeHeader(file);
	writeMacros(file);
}

void PostscriptWriter::write(Partition& clustering, std::string path) {
	TRACE("start ps write clustering");

	std::ofstream file;
	init(path, file);

	file << "0.00 0.00 0.00 c\n";   // 0.25 0.25 0.25 c 1.0 w

	writeClustering(clustering, file);

	if (! wrapAround) {
		file << "grestore\n";
	}
	file.close();
}

void PostscriptWriter::write(std::string path) {
	TRACE("start ps write");
	ClusteringGenerator gen;
	Partition allNone = gen.makeOneClustering(g);
	write(allNone, path);
}


} /* namespace NetworKit */
