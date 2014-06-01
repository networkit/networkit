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

	std::vector<count> graphSizes = {10000, 25000, 50000, 100000};
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
		printf("G(n = %lu, m = %lu), generation took %.1f s, memory used before shrink: %.1f KB, memory used after shrink: %.1f KB, relative saving: %.4f, shrinking took %lu %s\n",
			G.numberOfNodes(),
			G.numberOfEdges(),
			timerGen.elapsedMilliseconds() / 1000.0,
			before  / 1024.0,
			after / 1024.0,
			1.0 - (double) after / before,
			useMS ? timerShrink.elapsedMilliseconds() : timerShrink.elapsedMicroseconds(),
			useMS ? "ms" : "us");
	}
}

TEST_F(GraphMemoryBenchmark, shrinkToFitForGraphReader) {
	Aux::Timer timerRead, timerShrink;
	METISGraphReader reader;

	std::vector<std::string> graphFiles = {"../graphs/in-2004.graph", "../graphs/uk-2002.graph", "../graphs/uk-2007-05.graph"};
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
		printf("G(n = %lu, m = %lu), reading %s took %.1f ms, memory used before shrink: %.1f KB, memory used after shrink: %.1f KB, relative saving: %.4f, shrinking took %lu %s\n",
			G.numberOfNodes(),
			G.numberOfEdges(),
			graphFile.c_str(),
			timerRead.elapsedMilliseconds() / 1000.0,
			before / 1024.0,
			after / 1024.0,
			1.0 - (double) after / before,
			useMS ? timerShrink.elapsedMilliseconds() : timerShrink.elapsedMicroseconds(),
			useMS ? "ms" : "us");
	}
}

} /* namespace NetworKit */

#endif /*NOGTEST */
