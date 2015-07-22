/*
 * QuadTreePolarEuclidGTest.h
 *
 *  Created on: 22.07.2015
 *      Author: moritzl
 */

#ifndef QUADTREEPOLAREUCLIDGTEST_H_
#define QUADTREEPOLAREUCLIDGTEST_H_

#include <gtest/gtest.h>


namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity

class QuadTreePolarEuclidGTest : public testing::Test{
public:
	QuadTreePolarEuclidGTest();
	virtual ~QuadTreePolarEuclidGTest();
};

} /* namespace NetworKit */
#endif /* QUADTREEPOLAREUCLIDGTEST_H_ */
