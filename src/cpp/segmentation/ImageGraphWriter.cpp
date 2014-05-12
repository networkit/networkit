/*
 * ImageGraphWriter.cpp
 *
 *  Created on: 26.04.2014
 *      Author: Michael
 */

#include "ImageGraphWriter.h"

namespace NetworKit {

ImageGraphWriter::ImageGraphWriter() {
}

ImageGraphWriter::~ImageGraphWriter() {
}

void ImageGraphWriter::writePartitionsToSingleImage(Graph &graph, const Partition &partition, const std::string &pathToSourceImage) const {
	png::image<png::rgb_pixel> sourceImage(pathToSourceImage);
	png::image<png::rgb_pixel> partitionImage(sourceImage.get_width(), sourceImage.get_height());

	int counter = 0;
	for (auto subset : partition.subsetSizeMap()) {
		std::set<index> partitionMembers = partition.getMembers(subset.first);

		// calculate mean color
		std::vector<int> color(3, 0);
		for (auto node : partitionMembers) {
			Point<float> pixel = graph.getCoordinate(node);
			png::rgb_pixel pixelColor = sourceImage.get_pixel(pixel[0], pixel[1]);
			color[0] += pixelColor.red;
			color[1] += pixelColor.green;
			color[2] += pixelColor.blue;
		}

		color[0] = color[0] / partitionMembers.size();
		color[1] = color[1] / partitionMembers.size();
		color[2] = color[2] / partitionMembers.size();


		png::rgb_pixel pixelColor(color[0], color[1], color[2]);
		for (auto node : partitionMembers) {
			Point<float> pixel = graph.getCoordinate(node);
			partitionImage.set_pixel(pixel[0], pixel[1], pixelColor);
		}

		std::string path = pathToSourceImage.substr(0, pathToSourceImage.find_last_of(".")) + "_partitioned.png";
		partitionImage.write(path);

		counter++;
	}
}

void ImageGraphWriter::writePartitions(Graph &graph, const Partition &partition, const std::string &pathToSourceImage) const {
	png::image<png::rgb_pixel> sourceImage(pathToSourceImage);

	int counter = 0;
	for (auto subset : partition.subsetSizeMap()) {
		png::image<png::rgb_pixel> partitionImage = createImageForPartition(graph, partition.getMembers(subset.first), sourceImage);

		std::stringstream stream;
		stream << "_partition_";
		stream << counter;
		stream << ".png";

		std::string path = pathToSourceImage.substr(0, pathToSourceImage.find_last_of(".")) + stream.str();
		partitionImage.write(path);

		counter++;
	}
}

png::image<png::rgb_pixel> ImageGraphWriter::createImageForPartition(Graph &graph, const std::set<index> partitionMembers, const png::image<png::rgb_pixel> &sourceImage) const {
	png::image<png::rgb_pixel> image(sourceImage.get_width(), sourceImage.get_height());
	for (auto node : partitionMembers) {
		Point<float> pixel = graph.getCoordinate(node);
		image.set_pixel(pixel[0], pixel[1], sourceImage.get_pixel(pixel[0], pixel[1]));
	}

	return image;
}

} /* namespace NetworKit */
