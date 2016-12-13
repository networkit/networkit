/*
 * RasterReader.h
 *
 *  Created on: 20.07.2015
 *      Author: Henning
 */

#ifndef RASTERREADER_H_
#define RASTERREADER_H_

#include <utility>
#include <vector>
#include <string>

namespace NetworKit {

/**
 * Reader for NASA raster data of population.
 * Coordinate transformation into the unit square takes place.
 * Within each grid cell, points get a random location.
 */
class RasterReader {
private:
	double normalizationFactor;

public:
	/**
	 * @param[in] Factor by which population per grid cell should be multiplied
	 */
	RasterReader(double normalizationFactor);
	virtual ~RasterReader() = default;

	/**
	 * @param[in] path Path and name of file to be read.
	 * @return Vectors of x- and y-coordinates of the points.
	 */
	std::pair<std::vector<double>, std::vector<double> > read(const std::string& path);
};

} /* namespace NetworKit */
#endif /* RASTERREADER_H_ */
