/*
 * LAMGGTest.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael
 */

#ifndef LAMGGTEST_H_
#define LAMGGTEST_H_

#include "gtest/gtest.h"

#include "../../algebraic/Vector.h"
#include "../../io/METISGraphReader.h"
#include "../../io/METISGraphWriter.h"
#include "../../generators/BarabasiAlbertGenerator.h"
#include "../../components/ConnectedComponents.h"
#include "../../structures/Partition.h"

using namespace std;

namespace NetworKit {

class LAMGGTest : public testing::Test {
protected:
	const vector<string> GRAPH_INSTANCES = {"input/jazz.graph", "input/power.graph", "input/wing.graph"};

	Vector randZeroSum(const Graph& graph, size_t seed) const;
	Vector randVector(count dimension, double lower, double upper) const;
};

} /* namespace NetworKit */

#endif /* LAMGGTEST_H_ */
