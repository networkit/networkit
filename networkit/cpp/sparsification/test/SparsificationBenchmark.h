/*
 * SparsificationBenchmark.h
 *
 *  Created on: 31.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#ifndef SparsificationBENCHMARK_H_
#define SparsificationBENCHMARK_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../auxiliary/Timer.h"


namespace NetworKit {

class SparsificationBenchmark: public testing::Test {
protected:
	int64_t n;
public:
	SparsificationBenchmark();
	virtual ~SparsificationBenchmark();
	Graph makeCompleteGraph(count n);
};

} /* namespace NetworKit */
#endif /* SparsificationBENCHMARK_H_ */

#endif /*NOGTEST */
