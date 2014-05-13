/*
 * GraphMemoryBenchmark.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "GraphMemoryBenchmark.h"

#include "../Graph.h"
#include "../../auxiliary/Timer.h"
#include "../GraphGenerator.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

TEST_F(GraphMemoryBenchmark, shrinkToFitForErdosRenyi) {
	GraphGenerator gen;
	Aux::Timer timerGen, timerShrink;
	
	std::vector<count> graphSizes = {1000, 5000, 10000};
	for (auto& n : graphSizes) {
		timerGen.start();
		Graph G = gen.makeErdosRenyiGraph(n, 0.01);
		timerGen.stop();

		count before = G.getMemoryUsage();
		timerShrink.start();
		G.shrinkToFit();
		timerShrink.stop();
		count after = G.getMemoryUsage();

		printf("G(n = %lu, m = %lu), generation took %lu ms, memory used before shrink: %lu, memory used after shrink: %lu, relative saving: %f, shrinking took %lu us\n",
			G.numberOfNodes(),
			G.numberOfEdges(),
			timerGen.elapsedMilliseconds(),
			before,
			after,
			1.0 - (double) after / before,
			timerShrink.elapsedMicroseconds());
	}
}

TEST_F(GraphMemoryBenchmark, shrinkToFitForGraphReader) {
	Aux::Timer timerRead, timerShrink;
	METISGraphReader reader;

	timerRead.start();
	Graph G = reader.read("input/caidaRouterLevel.graph");
	timerRead.stop();

	count before = G.getMemoryUsage();
	timerShrink.start();
	G.shrinkToFit();
	timerShrink.stop();
	count after = G.getMemoryUsage();

	printf("G(n = %lu, m = %lu), reading took %lu ms, memory used before shrink: %lu, memory used after shrink: %lu, relative saving: %f, shrinking took %lu us\n",
		G.numberOfNodes(),
		G.numberOfEdges(),
		timerRead.elapsedMilliseconds(),
		before,
		after,
		1.0 - (double) after / before,
		timerShrink.elapsedMicroseconds());
}


} /* namespace NetworKit */

#endif /*NOGTEST */
