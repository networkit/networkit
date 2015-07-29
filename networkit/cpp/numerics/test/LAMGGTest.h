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
#include "../../properties/ConnectedComponents.h"
#include "../../structures/Partition.h"

using namespace std;

namespace NetworKit {

class LAMGGTest : public testing::Test {
protected:
	const vector<string> GRAPH_INSTANCES = {"walshaw/wing.graph"};

	const std::vector<count> grid2DSizes = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
	const std::vector<count> grid3DSizes = {2, 4, 8, 16, 32, 64, 128};
	const std::vector<count> barabasiSizes = {25, 100, 500, 1000, 5000, 10000, 50000};

	Vector randZeroSum(const Graph &graph, size_t seed) const;
	Vector randVector(count dimension, double lower, double upper) const;
	Vector capacitanceProblem(const Graph &graph) const;

	void genGrid2D(count n) const;
	void genGrid3D(count n) const;
	void genBarabasi(count n, count attachment = 4) const;

};

} /* namespace NetworKit */

#endif /* LAMGGTEST_H_ */
