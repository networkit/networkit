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
#include "../../auxiliary/Log.h"

#include "../HyperbolicGenerator.h"
#include "../DynamicHyperbolicGenerator.h"
#include "../../graph/GraphGenerator.h"
#include "../BarabasiAlbertGenerator.h"
#include "../../graph/GraphBuilder.h"

namespace NetworKit {


TEST_F(GeneratorsBenchmark, benchmarkGraphBuilder) {
	// parameters for Erdös-Renyi
	count n = 25000;
	double p = 0.001;
	count m_expected = p * n * (n + 1) / 2;

	Graph G;
	GraphBuilder builder;

	// prepare a random generator for each possible thread
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
	
	// half parallel way
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
		G = builder.toGraph(false);
	});
	m_actual = G.numberOfEdges();
	EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	std::cout << "parallelForNodePairs + toGraphSequentiel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";
	// printf("parallelForNodePairs + toGraphSequentiel:\t\t%" PRIu64 " + %" PRIu64 " = %" PRIu64 " ms\n", t1, t2, t1 + t2);

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
	std::cout << "parallelForNodePairs + toGraphParallel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";

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
	// std::cout << "forNodePairs + Graph.addEdge:\t\t\t\t" << t1 << " ms\n";
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

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGenerator) {
	count n = 500000;
	HyperbolicGenerator gen(n,1,1);
	Graph G = gen.generate();
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkDynamicHyperbolicGenerator) {
	count n = 10000;
	count nSteps = 100;
	//(count n, double initialFactor = 1, double alpha = 1, double stretch = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator dyngen(n,0,1,1,0,0.01,0);
	dyngen.generate(nSteps);


}

} /* namespace NetworKit */

#endif /*NOGTEST */
