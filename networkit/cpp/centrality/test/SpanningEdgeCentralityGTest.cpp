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

	Spanning sp(G);

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

	std::string graphFiles[2] = {"input/jazz.graph", "input/power.graph"};

	for (auto graphFile: graphFiles) {
		Graph G = reader.read(graphFile);
		G.indexEdges();
		Aux::Timer timer;
		Spanning exact(G);
		Spanning cen(G);

		timer.start();
		exact.run();
		timer.stop();
		INFO("spanning edge centrality ranking time: ", timer.elapsedTag());

		timer.start();
		cen.runApproximation();
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



		timer.start();
		cen.runTreeApproximation();
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
		cen.runPseudoTreeApproximation();
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

TEST_F(SpanningEdgeCentralityGTest, benchSpanning) {
	METISGraphReader reader;
	string benchFolder = "";
	std::ofstream output("/Users/Michael/Downloads/SNAP/spanningBenchmark" + "_0.1");
	if (!output.is_open()) {
		ERROR("Could not open output file");
		return;
	}

	Aux::Timer t;
	for (const string &instance : instances) {
		Graph G = reader.read(instance);
		G.indexEdges();


		Spanning sp(G);
		t.start();
		sp.runApproximation();
		t.stop();

		output << instance << "\t" << sp.getSetupTime() << "\t" << t.elapsedMilliseconds() << std::endl;
		output.flush();

		if (instance.find_first_of("CA-GrQc") != instance.npos) {
			std::ofstream scoreOutput(benchFolder + "_GrQC_scores_LAMG");
			if (scoreOutput.is_open()) {
				scoreOutput << std::setprecision(16);
				G.forEdges([&](node u, node v, edgeid e) {
					scoreOutput << u << " " << v << " " << sp.score(e) << std::endl;
				});
				scoreOutput.close();
			}
		}
	}

	output.close();
}

} /* namespace NetworKit */
