/*
 * BasicsBenchmark.h
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef BASICSBENCHMARK_H_
#define BASICSBENCHMARK_H_

#include <gtest/gtest.h>
#include <vector>

#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

class BasicsBenchmark: public testing::Test {
public:
	BasicsBenchmark();
	virtual ~BasicsBenchmark();
};

} /* namespace NetworKit */
#endif /* BASICSBENCHMARK_H_ */
