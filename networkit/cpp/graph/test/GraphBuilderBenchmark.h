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
#include <functional>
#include "../../auxiliary/Timer.h"

namespace NetworKit {

class GTestBenchmark : public testing::Test {
public:
	static void measureInMs(std::function<int()> func, int iterations = 1) {
		Aux::Timer timer;
		std::vector<unsigned long long> runningTimes(iterations);
		for (int i = 0; i < iterations; i++) {
			printf("Iteration %d of %d ...", i + 1, iterations);
			
			timer.start();
			int x = func();
			timer.stop();

			runningTimes[i] = timer.elapsedMilliseconds();

			printf("done (took %llu ms, result: %d).\n", runningTimes[i], x);
		}

		long long unsigned sum = 0;
		long long unsigned min = runningTimes[0];
		long long unsigned max = runningTimes[0];
		for (auto t : runningTimes) {
			sum += t;
			if (t < min) min = t;
			if (t > max) max = t;
		}

		printf("Iterations: %d, average runtime: %1.1f ms, fastest run: %llu ms, slowest run: %llu ms\n", iterations, (double) sum / iterations, min, max);
	}

	template <typename L>
	uint64_t timeOnce(L f) {
		// TODO should be moved somewhere else (Benchmark parent class or the Timer class itself)
		Aux::Timer timer;
		timer.start();
		f();
		timer.stop();
		return timer.elapsedMilliseconds();
	}
};

class GraphBuilderBenchmark : public GTestBenchmark {
public:
	GraphBuilderBenchmark() = default;
};

} /* namespace NetworKit */

#endif /* GRAPH_BUILDER_BENCHMARK_H_ */

#endif /*NOGTEST */
