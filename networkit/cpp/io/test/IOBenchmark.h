/*
 * IOBenchmark.h
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef IOBENCHMARK_H_
#define IOBENCHMARK_H_

#include <gtest/gtest.h>
#include <vector>
#include <string>

#include "../../auxiliary/Log.h"
#include "../../auxiliary/Timer.h"
#include "../METISGraphReader.h"

using std::vector;
using std::string;

namespace NetworKit {

class IOBenchmark: public testing::Test {
public:
	IOBenchmark() = default;
	virtual ~IOBenchmark() = default;

	static void convertToHeatMap(vector<bool> &infected, vector<double> &xcoords, vector<double> &ycoords, string filename, double resolution=1);
};

} /* namespace NetworKit */
#endif /* IOBENCHMARK_H_ */

#endif /* NOGTEST */
