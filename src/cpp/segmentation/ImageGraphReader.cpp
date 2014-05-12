/*
 * ImageGraphReader.cpp
 *
 *  Created on: 24.04.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "ImageGraphReader.h"
#include "../auxiliary/Log.h"
#include <time.h>

namespace NetworKit {

ImageGraphReader::ImageGraphReader() {
}

ImageGraphReader::~ImageGraphReader() {
}

Graph ImageGraphReader::read(const std::string &path, const double sigmaI, const double sigmaX, const int radius) const {
	std::ifstream imageStream(path, std::ios::binary);
	if (!imageStream.is_open()) {
		ERROR("invalid image file: ", path);
		throw std::runtime_error("invalid image file");
	}

	// check color type of image
	png::reader<std::istream> reader(imageStream);
	reader.read_info();
	png::color_type colorType = reader.get_color_type();

	// load image
	png::image<png::rgb_pixel> image(path);

	int width = image.get_width();
	int height = image.get_height();

	// initialize graph
	Graph graph(width * height, true);
	graph.initCoordinates();

	// add edges
	int rasterRad = radius - 1;
	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < height; ++j) {
			graph.setCoordinate(j * width + i, Point<float>(i, j));
			for (int v = std::max(0, (int) j - rasterRad); v <= j; ++v) {
				for (int h = v == j? i : std::max(0, i - rasterRad); h <= std::min(i + rasterRad, (int) width - 1); ++h) {
					std::vector<int> v1 = {i, j};
					std::vector<int> v2 = {h, v};
					graph.addEdge(j * width + i, v * width + h, getWeight(image, v1, v2, colorType, sigmaI, sigmaX));
				}
			}
		}
	}

	return graph;
}


double ImageGraphReader::inverseGammasRGB(const int &ic) const {
	double c = ic / 255.0;
	if ( c <= 0.04045 )
		return c / 12.92;
	else
		return std::pow(((c + 0.055) / (1.055)), 2.4);
}

int ImageGraphReader::gammasRGB(double v) const {
	if(v <= 0.0031308)
		v *= 12.92;
	else
		v = 1.055 * std::pow(v, 1.0 / 2.4) - 0.055;
	return int(v * 255 + 0.5);
}

int ImageGraphReader::gray(int r, int g, int b) const {
	return gammasRGB(rY * inverseGammasRGB(r) + gY * inverseGammasRGB(g) + bY * inverseGammasRGB(b));
}

double ImageGraphReader::getWeight(const png::image<png::rgb_pixel> &image, std::vector<int> &v1, std::vector<int> &v2, const png::color_type &colorType, const double sigmaI, const double sigmaX) const {
	png::rgb_pixel p1 = image.get_pixel(v1[0], v1[1]);
	png::rgb_pixel p2 = image.get_pixel(v2[0], v2[1]);
	double feature;
	if (colorType == png::color_type_gray || colorType == png::color_type_gray_alpha) {
		feature = std::pow(std::abs(gray(p1.red, p1.green, p1.blue) - gray(p2.red, p2.green, p2.blue)), 2);
	} else {
		std::vector<int> f1 = {p1.red, p1.green, p1.blue};
		std::vector<int> f2 = {p2.red, p2.green, p2.blue};

		// calculate f1 - f2 (and store it in f1)
		f1[0] = f1[0] - f2[0];
		f1[1] = f1[1] - f2[1];
		f1[2] = f1[2] - f2[2];

		// calculate scalar product of f1 (which is (f1 - f2))
		feature = f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2];
	}

	// calculate v1 - v2 (and store it in v1)
	v1[0] = v1[0] - v2[0];
	v1[1] = v1[1] - v2[1];

	// calculate scalar product of v1 (which is (v1 - v2))
	double sqDist = v1[0] * v1[0] + v2[0] * v2[0];

	return std::pow(M_E, -feature / (sigmaI * sigmaI)) * std::pow(M_E, -sqDist / (sigmaX * sigmaX));
}


} /** namespace NetworKit **/

