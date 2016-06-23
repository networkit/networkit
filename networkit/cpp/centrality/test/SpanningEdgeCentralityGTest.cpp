/*
 * SpanningEdgeCentralityGTest.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: Michael
 */

#include "SpanningEdgeCentralityGTest.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

#include <fstream>
#include <iomanip>


namespace NetworKit {

TEST_F(SpanningEdgeCentralityGTest, testOnToyGraph) {
	/* Graph:
		    0    3
		     \  / \
		      2    5
		     /  \ /
		    1    4
	 */
	count n = 6;
	Graph G(n, false, false);
	G.indexEdges();


	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);

	SpanningEdgeCentrality sp(G);

	sp.run();
	EXPECT_NEAR(1.0, sp.score(0), 1e-5);
	EXPECT_NEAR(1.0, sp.score(1), 1e-5);
	EXPECT_NEAR(0.75, sp.score(2), 1e-5);
	EXPECT_NEAR(0.75, sp.score(3), 1e-5);
	EXPECT_NEAR(0.75, sp.score(4), 1e-5);
	EXPECT_NEAR(0.75, sp.score(5), 1e-5);
}

TEST_F(SpanningEdgeCentralityGTest, testSpanningOnSmallGraphs) {
	METISGraphReader reader;

	std::string graphFiles[2] = {"input/karate.graph", "input/tiny_01.graph"};

	for (auto graphFile: graphFiles) {
		Graph G = reader.read(graphFile);
		G.indexEdges();
		Aux::Timer timer;
		SpanningEdgeCentrality exact(G);
		SpanningEdgeCentrality cen(G);

		timer.start();
		exact.run();
		timer.stop();
		INFO("spanning edge centrality ranking time: ", timer.elapsedTag());

		timer.start();
		cen.runParallelApproximation();
		timer.stop();
		INFO("approx spanning edge centrality time: ", timer.elapsedTag());

		double error = 0.0;
		G.forEdges([&](node u, node v, edgeid e) {
			double relError = fabs(cen.score(e) - exact.score(e));
			if (fabs(exact.score(e)) > 1e-9) relError /= exact.score(e);
			error += relError;
		});
		error /= G.numberOfEdges();
		INFO("Avg. relative error: ", error);

		count reps = 500;

		timer.start();
		cen.runTreeApproximation(reps);
		timer.stop();
		INFO("tree approx spanning edge centrality time: ", timer.elapsedTag());

		error = 0.0;
		G.forEdges([&](node u, node v, edgeid e) {
			double relError = fabs(cen.score(e) - exact.score(e));
			if (fabs(exact.score(e)) > 1e-9) relError /= exact.score(e);
			error += relError;
		});
		error /= G.numberOfEdges();
		INFO("Avg. relative error: ", error);


		timer.start();
		cen.runPseudoTreeApproximation(reps);
		timer.stop();
		INFO("pseudo tree approx spanning edge centrality time: ", timer.elapsedTag());

		error = 0.0;
		G.forEdges([&](node u, node v, edgeid e) {
			double relError = fabs(cen.score(e) - exact.score(e));
			if (fabs(exact.score(e)) > 1e-9) relError /= exact.score(e);
			error += relError;
		});
		error /= G.numberOfEdges();
		INFO("Avg. relative error: ", error);
	}
}

} /* namespace NetworKit */
