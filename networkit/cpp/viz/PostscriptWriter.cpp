/*
 * PostscriptWriter.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include "PostscriptWriter.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

PostscriptWriter::PostscriptWriter(bool isTorus) : wrapAround(isTorus) {
	numColors = 24;

	// set colors in RGB format, where 0.0 is no color and 1.0 is full color
	// (1.0 would mean 255 in a 3x8 bit color scheme)
	psColor = { {1.0, 0.0, 0.0},
		{	1.0, 0.5, 0.0}, {1.0, 1.0, 0.0}, {0.5, 1.0, 0.0}, {0.0, 1.0,
			0.0}, {0.0, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.5, 1.0},
		{	0.0, 0.0, 1.0}, {0.5, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 0.0,
			0.5},

		{	0.6, 0.0, 0.0}, {0.6, 0.3, 0.0}, {0.6, 0.6, 0.0}, {0.3, 0.6,
			0.0}, {0.0, 0.6, 0.0}, {0.0, 0.6, 0.3}, {0.0, 0.6, 0.6},
		{	0.0, 0.3, 0.6}, {0.0, 0.0, 0.6}, {0.3, 0.0, 0.6}, {0.6, 0.0,
			0.6}, {0.6, 0.0, 0.3}};

	// bounding box size
	ps_size = {1020.0, 1020.0};
}

void PostscriptWriter::writeHeader(std::ofstream& file) {
	/* Header */
	if (wrapAround) {
		file << "%!PS-Adobe-3.0 EPSF-3.0\n";
	} else {
		file << "%!PS-Adobe-1.0\n";
	}
	file << "%%Title: NetworKit visualization\n";
	file << "%%BoundingBox: 0.000 0.000 " << ps_size[0] << " " << ps_size[1] << "\n";
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
void PostscriptWriter::writeClustering(Graph& g, Partition& clustering, std::ofstream& file)
{
	/////////////////////////////////
	// bounding box adjustment
	/////////////////////////////////
	ps_min = {g.minCoordinate(0), g.minCoordinate(1)};
	ps_max = {g.maxCoordinate(0), g.maxCoordinate(1)};
	Point<float> ps_stretch = {ps_size[0] - 2 * ps_border[0], ps_size[1] - 2 * ps_border[1]};

	TRACE("min: ", ps_min.toCsvString());
	TRACE("max: ", ps_max.toCsvString());
	TRACE("stretch: ", ps_stretch.toCsvString());

	auto adjustToBoundingBox([&](Point<float> p) {
		for (index c = 0; c < 2; ++c) {
			p[c] -= ps_min[c];
			p[c] *= ps_stretch[c] / (ps_max[c] - ps_min[c]);
			p[c] += ps_border[c];
		}
//		TRACE("New coordinate: ", p.toCsvString());
		return p;
	});


	/////////////////////////////////
	// wrap-around adjustment
	/////////////////////////////////
	auto adjust1([&](float& val) { // TODO: externalize constants
		if (val > 500.0f) {
			val -= 1000.0f;
		}
		else if (val < -500.0f) {
			val += 1000.0f;
		}
	});

	auto adjustWrapAround([&](Point<float>& diff) {
		adjust1(diff[0]);
		adjust1(diff[1]);
	});


	// draw edges
	TRACE("start edge loop in writeClustering, wrapAround? ", wrapAround);
	TRACE("num edges: ", g.numberOfEdges());

	g.forEdges([&](node u, node v) {

		// set edge color
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

		// set edge start and end point
		Point<float> start = adjustToBoundingBox(g.getCoordinate(u));
		Point<float> end = adjustToBoundingBox(g.getCoordinate(v));
		Point<float> diff = {end[0] - start[0], end[1] - start[1]};
		if (wrapAround) {
			adjustWrapAround(diff);
		}
		end = start;
		end += diff;

		// write edge to file
		file << "p " << start.toSsvString() << " m " << end.toSsvString() << " l s\n";
	});


	// draw vertices
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

		Point<float> point = adjustToBoundingBox(g.getCoordinate(u));

		file << "p " << point.toSsvString() << " " << dotsize << " 0.00 360.00 a s\n";
//		TRACE("write coordinate to file: ", point[0], ", ", point[1]);
	});
}

void PostscriptWriter::init(std::string path, std::ofstream& file) {
	TRACE("start ps init");

	file.open(path.c_str());
	file.precision(3);
	file << std::fixed;

	writeHeader(file);
	writeMacros(file);
	file << "0.000 0.000 0.000 c\n";
}

void PostscriptWriter::write(Graph& g, Partition& clustering, std::string path) {
	TRACE("start ps write clustering");
	assert(g.getCoordinate(0).getDimensions() == 2);

	std::ofstream file;
	init(path, file);

	writeClustering(g, clustering, file);

	if (! wrapAround) {
		file << "grestore\n";
	}
	file.close();
}

void PostscriptWriter::write(Graph& g, std::string path) {
	TRACE("start ps write");
	ClusteringGenerator gen;
	Partition allNone = gen.makeOneClustering(g);
	write(g, allNone, path);
}


} /* namespace NetworKit */
