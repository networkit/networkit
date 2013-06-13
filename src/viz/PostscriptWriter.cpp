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
void PostscriptWriter::writeClustering(Clustering& clustering,
		std::ofstream& file) {
	auto adjust([&](float& val) {
		val += 10.0;
	});

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

		float startx = g.getCoordinate(u, 0);
		float starty = g.getCoordinate(u, 1);
		float endx = g.getCoordinate(v, 0);
		float endy = g.getCoordinate(v, 1);

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
	float dotsize = 3.6;
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

		float x = g.getCoordinate(u, 0);
		float y = g.getCoordinate(u, 1);

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

//		DEBUG("x: " << x);
//		DEBUG("y: " << y);

	});
}

void PostscriptWriter::init(std::string path, std::ofstream& file) {
	file.open(path.c_str());
	if (true || wrapAround) { // FIXME
		file.precision(3);
		// adjust coordinates for postscript output
		g.forNodes([&](node u) {
			g.setCoordinate(u, 0, 1000.0 * g.getCoordinate(u, 0));
			g.setCoordinate(u, 1, 1000.0 * g.getCoordinate(u, 1));
		});
	} else {
		file.precision(2);
	}
	file << std::fixed;
}

void PostscriptWriter::write(Clustering& clustering, std::string path) {
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
	std::ofstream file;
	init(path, file);
	writeHeader(file);
	writeMacros(file);

	ClusteringGenerator gen;
	Clustering allNone = gen.makeOneClustering(g);
	writeClustering(allNone, file);
	if (!wrapAround) {
		file << "grestore\n";
	}

	file.close();
}

void PostscriptWriter::writeAlgebraicDistances(std::string path, const AlgebraicDistances& ad) {
//	int minvw = INT_MAX, maxvw = INT_MIN;
//	int minew = INT_MAX, maxew = INT_MIN;
	float minload = (float) INT_MAX, maxload = (float) INT_MIN;
	float mindist = (float) INT_MAX, maxdist = (float) INT_MIN;
	float dotsize = 3.0;
	float red, green, blue, range;

	FILE * ps = NULL;

	printf("  writing postscript file '%s'", path.c_str());

	if ((ps = fopen(path.c_str(), "w")) == NULL) {
		throw std::runtime_error("cannot open file");
		return;
	}

	/* Header */
	fprintf(ps, "%%!PS-Adobe-1.0\n");
	fprintf(ps, "%%%%Title: %s\n", path.c_str());
	fprintf(ps, "%%%%BoundingBox: 0 0 %i, %i\n", (int) ps_sizex,
			(int) ceil(ps_sizey) + 40);
	fprintf(ps, "%%%%EndComments\n");
	fprintf(ps, "%%%%EndProlog\n");
	fprintf(ps, "gsave\n");

	/* Macros */
	fprintf(ps, "/p {newpath} bind def\n");
	fprintf(ps, "/m {moveto} bind def\n");
	fprintf(ps, "/r {rmoveto} bind def\n");
	fprintf(ps, "/l {lineto} bind def\n");
	fprintf(ps, "/n {rlineto} bind def\n");
	fprintf(ps, "/c {setrgbcolor} bind def\n");
	fprintf(ps, "/s {stroke} bind def\n");
	fprintf(ps, "/w {setlinewidth} bind def\n");
	fprintf(ps, "/h {show} bind def\n");
	fprintf(ps, "/a {arc closepath fill} bind def\n");
	fprintf(ps, "/b {closepath eofill} bind def\n");
	fprintf(ps, "0.00 0.00 0.00 c\n");

	minload = 0.0f;

	g.forNodes([&](node i) {
		// max and min load
		double load = ad.geometricMeanLoad(i);
		if (load > maxload) {
			maxload = load;
		}
		if (load < minload) {
			minload = load;
		}

		// max and min dist
		g.forNeighborsOf(i, [&](node j) {
			double dist = ad.algdist(i, j);
			if (dist > maxdist)
				maxdist = dist;
			if (dist < mindist)
				mindist = dist;
		});
	});

	/* Kanten zeichnen */
	g.forNodes([&](node i) {
		g.forNeighborsOf(i, [&](node j) {
			if (mindist == maxdist || ad.geometricMeanLoad(j) == 0.0)
				fprintf(ps, "0.50 0.50 0.50 c ");
			else {
				range = fabs(ad.algdist(i, j) - mindist) / (maxdist - mindist);
				if (range > 0.8) {
					red = 1.0;
					green = 1.0 - (range - 0.8) * 5.0;
					blue = 0.0;
				} else if (range > 0.6) {
					red = (range - 0.6) * 5.0;
					green = 1.0;
					blue = 0.0;
				} else if (range > 0.4) {
					red = 0.0;
					green = 1.0;
					blue = 1.0 - (range - 0.4) * 5.0;
				} else if (range > 0.2) {
					red = 0.0;
					green = (range - 0.2) * 5.0;
					blue = 1.0;
				} else {
					red = 0.0;
					green = 0.0;
					blue = range * 5.0;
				}
				fprintf(ps, "%.2f %.2f %.2f c ", red, green, blue);
			}

			fprintf(ps, "p %.2f %.2f m %.2f %.2f l s\n",
					(g.getCoordinate(i, 0) - ps_minx) * ps_scale + ps_borderx,
					(g.getCoordinate(i, 1) - ps_miny) * ps_scale + ps_bordery,
					(g.getCoordinate(j, 0) - ps_minx) * ps_scale + ps_borderx,
					(g.getCoordinate(j, 1) - ps_miny) * ps_scale + ps_bordery);
		});
	});

	/* Knoten zeichnen */
	g.forNodes([&](node i) {
		if (minload == maxload)
			fprintf(ps, "0.50 0.50 0.50 c ");
		else {
			range = (ad.geometricMeanLoad(i) - minload) / (maxload - minload);
			if (range > 0.8) {
				red = 1.0;
				green = 1.0 - (range - 0.8) * 5.0;
				blue = 0.0;
			} else if (range > 0.6) {
				red = (range - 0.6) * 5.0;
				green = 1.0;
				blue = 0.0;
			} else if (range > 0.4) {
				red = 0.0;
				green = 1.0;
				blue = 1.0 - (range - 0.4) * 5.0;
			} else if (range > 0.2) {
				red = 0.0;
				green = (range - 0.2) * 5.0;
				blue = 1.0;
			} else {
				red = 0.0;
				green = 0.0;
				blue = range * 5.0;
			}
			fprintf(ps, "%.2f %.2f %.2f c ", red, green, blue);
		}

		fprintf(ps, "p %.2f %.2f %.2f 0.00 360.00 a\n",
				(g.getCoordinate(i, 0) - ps_minx) * ps_scale + ps_borderx,
				(g.getCoordinate(i, 1) - ps_miny) * ps_scale + ps_bordery, dotsize);
	});

	fprintf(ps, "grestore\n");
	fclose(ps);
}

} /* namespace NetworKit */
