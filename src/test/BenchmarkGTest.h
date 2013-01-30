/*
 * BenchmarkGTest.h
 *
 *  Created on: 30.01.2013
 *      Author: cls
 */

#ifndef BENCHMARKGTEST_H_
#define BENCHMARKGTEST_H_

#include <gtest/gtest.h>
#include <vector>

#include "../aux/Log.h"
#include "../aux/Timer.h"

namespace EnsembleClustering {

class BenchmarkGTest: public testing::Test {
public:
	BenchmarkGTest();
	virtual ~BenchmarkGTest();
};

} /* namespace EnsembleClustering */
#endif /* BENCHMARKGTEST_H_ */
