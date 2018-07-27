/*
 * GeometricGTest.h
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#ifndef GEOMETRICGTEST_H_
#define GEOMETRICGTEST_H_

#include <gtest/gtest.h>

namespace NetworKit {

class GeometricGTest: public testing::Test {
public:
	GeometricGTest() = default;
	virtual ~GeometricGTest() = default;
};

} /* namespace NetworKit */
#endif /* GEOMETRICGTEST_H_ */
