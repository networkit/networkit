/*
 * GraphBuilderBenchmark.h
 *
 *  Created on: 04.12.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef GRAPH_BUILDER_BENCHMARK_H_
#define GRAPH_BUILDER_BENCHMARK_H_

#include <gtest/gtest.h>

#include "../GraphBuilder.h"
#include "../../auxiliary/Timer.h"

namespace NetworKit {

class GTestBenchmark : public testing::Test {
public:
	static void measureInMs(std::function<int()> func, int iterations = 1) {
		Aux::Timer timer;
		std::vector<uint64_t> runningTimes(iterations);
		for (int i = 0; i < iterations; i++) {
			printf("Iteration %d of %d ...", i + 1, iterations);
			
			timer.start();
			int x = func();
			timer.stop();

			runningTimes[i] = timer.elapsedMilliseconds();

			printf("done (took %lu ms, result: %d).\n", runningTimes[i], x);
		}

		uint64_t sum = 0, min = runningTimes[0], max = runningTimes[0];
		for (auto t : runningTimes) {
			sum += t;
			if (t < min) min = t;
			if (t > max) max = t;
		}

		printf("Iterations: %d, average runtime: %1.1f ms, fastest run: %lu ms, slowest run: %lu ms\n", iterations, (double) sum / iterations, min, max);
	}
};

class GraphBuilderBenchmark : public GTestBenchmark {
public:
	GraphBuilderBenchmark();
};

} /* namespace NetworKit */

#endif /* GRAPH_BUILDER_BENCHMARK_H_ */

#endif /*NOGTEST */