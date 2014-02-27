/*
 * PostscriptWriterGTest.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#ifndef POSTSCRIPTWRITERGTEST_H_
#define POSTSCRIPTWRITERGTEST_H_

#include <gtest/gtest.h>
#include <vector>

#include "../PostscriptWriter.h"
#include "../FruchtermanReingold.h"
#include "../MaxentStress.h"
#include "../MultilevelLayouter.h"
#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../io/ClusteringWriter.h"
#include "../../io/METISGraphReader.h"
#include "../../generators/PubWebGenerator.h"


namespace NetworKit {

class VizGTest : public testing::Test {
public:
	VizGTest();
	virtual ~VizGTest();
};

} /* namespace NetworKit */
#endif /* POSTSCRIPTWRITERGTEST_H_ */

#endif /*NOGTEST */
