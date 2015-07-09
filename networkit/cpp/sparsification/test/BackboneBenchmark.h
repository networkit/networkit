/*
 * BackboneBenchmark.h
 *
 *  Created on: 31.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#ifndef BACKBONEBENCHMARK_H_
#define BACKBONEBENCHMARK_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../auxiliary/Timer.h"


namespace NetworKit {

class BackboneBenchmark: public testing::Test {
protected:
	int64_t n;
public:
	BackboneBenchmark();
	virtual ~BackboneBenchmark();
	Graph makeCompleteGraph(count n);
};

} /* namespace NetworKit */
#endif /* BACKBONEBENCHMARK_H_ */

#endif /*NOGTEST */
