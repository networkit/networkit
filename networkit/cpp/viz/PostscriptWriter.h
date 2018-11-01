/*
 * PostscriptWriter.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef POSTSCRIPTWRITER_H_
#define POSTSCRIPTWRITER_H_

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <climits>

#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../community/ClusteringGenerator.h"

namespace NetworKit {

/**
 * @ingroup viz
 * EPS output of graphs with 2D coordinates
 */
class PostscriptWriter {

	struct float_triple {
		float r;
		float g;
		float b;
	};

protected:
	bool wrapAround;

	count numColors;
	std::vector<float_triple> psColor;

	Point<float> ps_size;
	Point<float> ps_border;
	Point<float> ps_min;
	Point<float> ps_max;

	void init(std::string filename, std::ofstream& file);
	void writeHeader(std::ofstream& file);
	void writeMacros(std::ofstream& file);
	void writeClustering(Graph& g, Partition& clustering, std::ofstream& file);

public:
	/**
	 * @param[in] isTorus Specifies whether the visualization square is treated as torus,
	 * i.e. with wrap-around boundaries (edge can leave the square and enter at the opposite
	 * side. By default, it is set to false.
	 */
	PostscriptWriter(bool isTorus = false);

	/**
	 * Outputs an EPS file with name @a filename of the graph @a g with 2D coordinates.
	 * The colors are chosen to visualize the specified @a clustering.
	 * @param[in] g Graph to be visualized.
	 * @param[in] clustering Clustering of the graph, visualized by different colors.
	 * @param[in] filename Name of file to write to.
	 */
	void write(Graph& g, Partition& clustering, std::string filename);

	/**
	 * Outputs an EPS file with name @a filename of the graph @a g with 2D coordinates.
	 * @param[in] g Graph to be visualized.
	 * @param[in] filename Name of file to write to.
	 */
	void write(Graph& g, std::string filename);
};

} /* namespace NetworKit */
#endif /* POSTSCRIPTWRITER_H_ */
