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

		bool useMS = timerShrink.elapsedMicroseconds() > 5000;
		printf("G(n = %lu, m = %lu), generation took %lu ms, memory used before shrink: %lu KB, memory used after shrink: %lu KB, relative saving: %f, shrinking took %lu %s\n",
			G.numberOfNodes(),
			G.numberOfEdges(),
			timerGen.elapsedMilliseconds(),
			before  / 1024,
			after / 1024,
			1.0 - (double) after / before,
			useMS ? timerShrink.elapsedMilliseconds() : timerShrink.elapsedMicroseconds(),
			useMS ? "ms" : "us");
	}
}

TEST_F(GraphMemoryBenchmark, shrinkToFitForGraphReader) {
	Aux::Timer timerRead, timerShrink;
	METISGraphReader reader;

	std::vector<std::string> graphFiles = {"input/astro-ph.graph", "input/caidaRouterLevel.graph", "../in-2004.graph"};
	for (auto& graphFile : graphFiles) {
		timerRead.start();
		Graph G = reader.read(graphFile);
		timerRead.stop();

		count before = G.getMemoryUsage();
		timerShrink.start();
		G.shrinkToFit();
		timerShrink.stop();
		count after = G.getMemoryUsage();

		bool useMS = timerShrink.elapsedMicroseconds() > 5000;
		printf("G(n = %lu, m = %lu), reading %s took %lu ms, memory used before shrink: %lu KB, memory used after shrink: %lu KB, relative saving: %f, shrinking took %lu %s\n",
			G.numberOfNodes(),
			G.numberOfEdges(),
			graphFile.c_str(),
			timerRead.elapsedMilliseconds(),
			before / 1024,
			after / 1024,
			1.0 - (double) after / before,
			useMS ? timerShrink.elapsedMilliseconds() : timerShrink.elapsedMicroseconds(),
			useMS ? "ms" : "us");
	}
}


} /* namespace NetworKit */

#endif /*NOGTEST */
