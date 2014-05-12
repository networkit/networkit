/*
 * ImageSegmentation.cpp
 *
 *  Created on: 27.04.2014
 *      Author: Michael
 */

#include "ImageSegmentation.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ImageSegmentation::ImageSegmentation() {
}

ImageSegmentation::~ImageSegmentation() {
}

void ImageSegmentation::createSegmentation(CommunityDetectionAlgorithm &algorithm, const std::string &path, bool allPartitionsToOneFile) {
	ImageGraphReader reader;
	ImageGraphWriter writer;

	TRACE("reading graph from: " , path);
	Graph graph = reader.read(path);

	TRACE("graph read, running clustering algorithm");
	Partition partition = algorithm.run(graph);

	TRACE("graph clustered, writing partitions to disk");
	if (allPartitionsToOneFile) {
		writer.writePartitionsToSingleImage(graph, partition, path);
	} else {
		writer.writePartitions(graph, partition, path);
	}
}


} /* namespace NetworKit */
