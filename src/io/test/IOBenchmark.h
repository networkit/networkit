/*
 * IOBenchmark.h
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef IOBENCHMARK_H_
#define IOBENCHMARK_H_

#include <gtest/gtest.h>

#include "../../aux/Log.h"
#include "../../aux/Timer.h"
#include "../METISGraphReader.h"

namespace EnsembleClustering {

class IOBenchmark: public testing::Test {
public:
	IOBenchmark();
	virtual ~IOBenchmark();
};

} /* namespace EnsembleClustering */
#endif /* IOBENCHMARK_H_ */
