/*
 * ImageGraphReaderGTest.h
 *
 *  Created on: 24.04.2014
 *      Author: Michael
 */

#ifndef IMAGEGRAPHREADERGTEST_H_
#define IMAGEGRAPHREADERGTEST_H_

#include "gtest/gtest.h"
#include "../ImageGraphReader.h"
#include "../../graph/Graph.h"

namespace NetworKit {

class ImageGraphReaderGTest : public testing::Test {
public:
	ImageGraphReaderGTest();
	virtual ~ImageGraphReaderGTest();
};

} /* namespace NetworKit */

#endif /* IMAGEGRAPHREADERGTEST_H_ */
