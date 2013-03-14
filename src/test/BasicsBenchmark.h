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

#include "../aux/Log.h"
#include "../aux/Timer.h"

namespace EnsembleClustering {

class BasicsBenchmark: public testing::Test {
public:
	BasicsBenchmark();
	virtual ~BasicsBenchmark();
};

} /* namespace EnsembleClustering */
#endif /* BASICSBENCHMARK_H_ */
