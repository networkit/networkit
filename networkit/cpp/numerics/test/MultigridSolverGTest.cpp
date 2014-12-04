/*
 * MultigridSolverGTest.cpp
 *
 *  Created on: 03.11.2014
 *      Author: Michael
 */

#include "MultigridSolverGTest.h"

using namespace std;

namespace NetworKit {

TEST_F(MultigridSolverGTest, tryOneCycle) {
	vector<Vector> rows;
	rows.push_back({10, -1, 2, 0});
	rows.push_back({-1, 11, -1, 3});
	rows.push_back({2, -1, 10, -1});
	rows.push_back({0, 3, -1, 8});
	Matrix A(rows);

	Vector b = {6, 25, -11, 15};
	Vector result = {0, 0, 0, 0};

	GaussSeidelRelaxation smoother;
	MatchingHierarchyBuilder hierarchyBuilder;
	MultigridSolver solver(1e-6, 3, 3, hierarchyBuilder, smoother);

	solver.setup(A);
	bool converged = solver.solve(b, result);

	EXPECT_TRUE(converged);
	EXPECT_EQ(1, std::round(result[0]));
	EXPECT_EQ(2, std::round(result[1]));
	EXPECT_EQ(-1, std::round(result[2]));
	EXPECT_EQ(1, std::round(result[3]));
}

TEST_F(MultigridSolverGTest, tryLarge) {
	METISGraphReader reader;
	GaussSeidelRelaxation smoother;
	MatchingHierarchyBuilder hierarchyBuilder;
	MultigridSolver solver(1e-3, 3, 3, hierarchyBuilder, smoother);

	for (const string &graph : GRAPH_INSTANCES) {
		Graph G = reader.read("input/" + graph);
		TRACE("setting up solver for graph ", graph);
		solver.setup(G);
		TRACE("DONE");

		Vector b = randZeroSum(G, 1234);
		Vector result(G.upperNodeIdBound(), 0.0);

		TRACE("Solving equation system");
		bool converged = solver.solve(b, result);
		TRACE("DONE");

		EXPECT_TRUE(converged);
	}
}

Vector MultigridSolverGTest::randZeroSum(const Graph& G, size_t seed) {
	mt19937 rand(seed);
	auto rand_value = uniform_real_distribution<double>(-1.0, 1.0);
	ConnectedComponents con(G);
	count n = G.numberOfNodes();
	con.run();
	Partition comps = con.getPartition();

	/* Fill each component randomly such that its sum is 0 */
	Vector b(n, 0.0);
	for (int id : comps.getSubsetIds()) {
		auto indexes = comps.getMembers(id);
		assert(!indexes.empty());
		double sum = 0.0;
		for (auto entry : indexes) {
			b[entry] = rand_value(rand);
			sum += b[entry];
		}
		b[*indexes.begin()] -= sum;
	}
	return b;
}

} /* namespace NetworKit */
