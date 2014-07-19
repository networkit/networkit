/*
 * GeneratorsBenchmark.cpp
 *
 *  Created on: May 29, 2013
 *      Author: forigem
 */

#ifndef NOGTEST

#include <omp.h>
#include <random>
#include <functional>

#include "GeneratorsBenchmark.h"
#include "../../graph/GraphGenerator.h"
#include "../BarabasiAlbertGenerator.h"
#include "../../graph/GraphBuilder.h"

namespace NetworKit {


TEST_F(GeneratorsBenchmark, benchmarkGraphBuilder) {
	// parameters for Erdös-Renyi
	count n = 100000;
	double p = 0.001;
	count m_expected = p * n * (n + 1) / 2;

	Graph G;
	GraphBuilder builder;

	// for (int powThread = 5; powThread < 10; powThread++) {
	// 	int threads = std::pow(2, powThread);
	// 	omp_set_num_threads(threads);
	// 	printf("\nthreads: %d\n", threads);

	// prepare a random generator for each possible thread
	// omp_set_num_threads(256);
	int maxThreads = omp_get_max_threads();
	std::vector< std::function<double()> > randomPerThread;
	std::random_device device;
	std::uniform_int_distribution<uint64_t> intDist;
	for (int tid = 0; tid < maxThreads; tid++) {
		auto seed = intDist(device);
		std::mt19937_64 gen(seed);
		std::uniform_real_distribution<double> dist{0.0, std::nexttoward(1.0, 2.0)};
		auto rdn = std::bind(dist, gen);
		randomPerThread.push_back(rdn);
	}

	count m_actual;
	uint64_t t1, t2;
	
	// // builder completely sequential
	// t1 = timeOnce([&]() {
	// 	builder = GraphBuilder(n);
	// 	builder.forNodePairs([&](node u, node v) {
	// 		if (randomPerThread[0]() <= p) {
	// 			builder.addEdge(u, v);
	// 		}
	// 	});
	// });
	// t2 = timeOnce([&]() {
	// 	G = builder.toGraph(false);
	// });
	// m_actual = G.numberOfEdges();
	// EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	// printf("forNodePairs + toGraphSequential:\t\t%lu + %lu = %lu ms\n", t1, t2, t1 + t2);

	// // parallel construction, but sequential toGraph
	// t1 = timeOnce([&]() {
	// 	builder = GraphBuilder(n);
	// 	builder.parallelForNodePairs([&](node u, node v) {
	// 		int tid = omp_get_thread_num();
	// 		double rdn = randomPerThread[tid]();
	// 		if (rdn <= p) {
	// 			builder.addEdge(u, v);
	// 		}
	// 	});
	// });
	// t2 = timeOnce([&]() {
	// 	G = builder.toGraph(false);
	// });
	// m_actual = G.numberOfEdges();
	// EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	// printf("parallelForNodePairs + toGraphSequential:\t%lu + %lu = %lu ms\n", t1, t2, t1 + t2);

	// fully parallel way
	m_actual = 0;
	t1 = timeOnce([&]() {
		builder = GraphBuilder(n);
		builder.parallelForNodePairs([&](node u, node v) {
			int tid = omp_get_thread_num();
			double rdn = randomPerThread[tid]();
			if (rdn <= p) {
				builder.addEdge(u, v);
			}
		});
	});
	t2 = timeOnce([&]() {
		G = builder.toGraph();
	});
	m_actual = G.numberOfEdges();
	EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	printf("parallelForNodePairs + toGraphParallel:\t\t%lu + %lu = %lu ms\n", t1, t2, t1 + t2);

	// old way
	// t1 = timeOnce([&]() {
	// 	G = Graph(n);
	// 	G.forNodePairs([&](node u, node v) {
	// 		if (randomPerThread[0]() <= p) {
	// 			G.addEdge(u, v);
	// 		}
	// 	});
	// });
	// m_actual = G.numberOfEdges();
	// EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	// printf("forNodePairs + Graph.addEdge:\t\t\t\t%lu ms\n", t1);

	// }
}

TEST_F(GeneratorsBenchmark, benchmarkBarabasiAlbertGenerator) {
	count k = 2;
	count nMax = 100000;
	count n0 = 2;

	BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());
}

} /* namespace NetworKit */

#endif /*NOGTEST */
