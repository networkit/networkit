/*
 * MultigridSolverGTest.h
 *
 *  Created on: 03.11.2014
 *      Author: Michael
 */

#ifndef NOGTEST

#ifndef MULTIGRIDSOLVERGTEST_H_
#define MULTIGRIDSOLVERGTEST_H_

#include "gtest/gtest.h"

#include "../../algebraic/Matrix.h"
#include "../../algebraic/Vector.h"
#include "../../io/METISGraphReader.h"
#include "../../properties/ConnectedComponents.h"
#include "../../structures/Partition.h"

#include "../MultigridSolver.h"
#include "../MatchingHierarchyBuilder.h"
#include "../GaussSeidelRelaxation.h"

using namespace std;

namespace NetworKit {

class MultigridSolverGTest : public testing::Test {
protected:
	const vector<string> GRAPH_INSTANCES = {"lesmis.graph", "power.graph", "celegans_metabolic.graph", "airfoil1.graph", "PGPgiantcompo.graph"};

	Vector randZeroSum(const Graph &graph, size_t seed);
};

} /* namespace NetworKit */

#endif /* MULTIGRIDSOLVERGTEST_H_ */

#endif
