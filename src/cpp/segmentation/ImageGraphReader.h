/*
 * ImageGraphReader.h
 *
 *  Created on: 24.04.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef IMAGEGRAPHREADER_H_
#define IMAGEGRAPHREADER_H_

#include "../algebraic/Vector.h"
#include "../graph/Graph.h"

#include "../imageio/png++/reader.hpp"
#include "../imageio/png++/image.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace NetworKit {

/**
 * The ImageGraphReader class reads png image files and creates Graph representations of them.
 *
 * The edge weight function is based on the paper "Normalized Cuts and Image Segmentation", Jianbo Shi
 * and Jitendra Malik, IEEE Transactions on Pattern Analysis and Machine Intelligence, vol.22, pp. 888-905,
 * 1997.
 */
class ImageGraphReader {
public:
	/** Default constructor. */
	ImageGraphReader();

	/** Destructor */
	virtual ~ImageGraphReader();

	/**
	 * Reads the png image located at @a path and constructs a Graph representation of it.
	 * @param path The file path to the png image.
	 * @param sigmaI Regulates the sensitivity to differences in brightness or color. Use higher values for less sensitivity.
	 * @param sigmaX Regulates the impact of more distant pixels. Use higher values for less impact.
	 * @param r Defines the radius of adjacent pixels to which edges are created. Edges are added to pixels with distance [0, r) from the center pixel.
	 * @return The constructed Graph from the image.
	 */
	Graph read(const std::string &path, const double sigmaI = 15, const double sigmaX = 75, const int radius = 4) const;

private:
	// sRGB luminance(Y) values
	const double rY = 0.212655;
	const double gY = 0.715158;
	const double bY = 0.072187;

	// Inverse of sRGB "gamma" function. (approx 2.2)
	double inverseGammasRGB(const int &ic) const;

	// sRGB "gamma" function (approx 2.2)
	int gammasRGB(double v) const;

	// gray value ("brightness")
	int gray(int r, int g, int b) const;

	// computes the weight of the edge between v1 and v2 according to the paper cited in the class description.
	double getWeight(const png::image<png::rgb_pixel> &image, std::vector<int> &v1, std::vector<int> &v2, const png::color_type &colorType, const double sigmaI, const double sigmaX) const;
};

} /** namespace NetworKit **/

#endif /* IMAGEGRAPHREADER_H_ */
