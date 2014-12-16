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
	count n = 100000;
	HyperbolicGenerator gen;
	Graph G = gen.generate(n,1,1);
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGeneratorWithSortedNodes) {
	count n = 100000;
	double s = 1.0;
	double alpha = 1.0;
	double t = 1.0;
	vector<double> angles(n);
	vector<double> radii(n);
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	//sample points randomly

	HyperbolicSpace::fillPoints(angles, radii, s, alpha);
	vector<index> permutation(n);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	//can probably be parallelized easily, but doesn't bring much benefit
	std::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> anglecopy(n);
	vector<double> radiicopy(n);

	#pragma omp parallel for
	for (index j = 0; j < n; j++) {
		anglecopy[j] = angles[permutation[j]];
		radiicopy[j] = radii[permutation[j]];
	}

	Graph G = HyperbolicGenerator().generate(anglecopy, radiicopy, r, R*t);
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGeneratorWithParallelQuadtree) {
	count n = 100000;
	double s = 1.0;
	Quadtree<index> quad(n,s);
	vector<double> angles;
	vector<double> radii;
	quad.trim();
	quad.sortPointsInLeaves();
	quad.reindex();
	quad.extractCoordinates(angles, radii);
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);

	HyperbolicGenerator gen;
	Graph G = gen.generate(angles, radii, quad, R);

	EXPECT_EQ(n, G.numberOfNodes());
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGeneratorWithSequentialQuadtree) {
	count n = 100000;
	double s = 1.0;
	double alpha = 1;

	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, s, alpha);
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	Quadtree<index> quad(r);

	for (index i = 0; i < n; i++) {
		quad.addContent(i, angles[i], radii[i]);
	}

	angles.clear();
	radii.clear();

	quad.trim();
	quad.sortPointsInLeaves();
	quad.reindex();
	quad.extractCoordinates(angles, radii);

	HyperbolicGenerator gen;
	Graph G = gen.generate(angles, radii, quad, R);

	EXPECT_EQ(n, G.numberOfNodes());
}

TEST_F(GeneratorsBenchmark, benchmarkDynamicHyperbolicGeneratorOnFactorGrowth) {
	count n = 10000;
	count nSteps = 100;
	//(count n, double initialFactor = 1, double alpha = 1, double stretch = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator dyngen(n,0,1,1,0,0.01,0);
	dyngen.generate(nSteps);
}

TEST_F(GeneratorsBenchmark, benchmarkDynamicHyperbolicGeneratorOnNodeMovement) {
	count n = 10000;
	count nSteps = 100;
	//(count n, double initialFactor = 1, double alpha = 1, double stretch = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator dyngen(n,1,1,1,0.5,0,0.02);
	dyngen.generate(nSteps);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
