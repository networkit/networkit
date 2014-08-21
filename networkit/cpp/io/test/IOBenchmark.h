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

#include "../../auxiliary/Log.h"
#include "../../auxiliary/Timer.h"
#include "../METISGraphReader.h"

namespace NetworKit {

class IOBenchmark: public testing::Test {
public:
	IOBenchmark();
	virtual ~IOBenchmark();
};

} /* namespace NetworKit */
#endif /* IOBENCHMARK_H_ */

#endif /* NOGTEST */
