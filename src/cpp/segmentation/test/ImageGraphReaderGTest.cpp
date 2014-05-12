/*
 * ImageGraphReaderGTest.cpp
 *
 *  Created on: 24.04.2014
 *      Author: Michael
 */

#include "ImageGraphReaderGTest.h"

namespace NetworKit {

ImageGraphReaderGTest::ImageGraphReaderGTest() {
}

ImageGraphReaderGTest::~ImageGraphReaderGTest() {
}

TEST(ImageGraphReaderGTest, tryReadImageTest) {
	ImageGraphReader reader;
	Graph graph = reader.read("input/segmentation/sample_1.png");

	ASSERT_EQ(16384u, graph.numberOfNodes());
	ASSERT_TRUE(graph.isWeighted());
}

} /* namespace NetworKit */
