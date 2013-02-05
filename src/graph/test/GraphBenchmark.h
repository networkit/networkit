/*
 * GraphBenchmark.h
 *
 *  Created on: 01.02.2013
 *      Author: cls
 */

#ifndef GRAPHBENCHMARK_H_
#define GRAPHBENCHMARK_H_

#include <gtest/gtest.h>

#include "../GraphGenerator.h"
#include "../../aux/Timer.h"

extern "C" {
#include "stinger.h"
}

namespace EnsembleClustering {

class GraphBenchmark: public testing::Test {
protected:
	int64_t n;
public:
	GraphBenchmark();
	virtual ~GraphBenchmark();
};

} /* namespace EnsembleClustering */
#endif /* GRAPHBENCHMARK_H_ */
