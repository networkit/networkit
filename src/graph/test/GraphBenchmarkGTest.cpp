/*
 * GraphBenchmarkGTest.cpp
 *
 *  Created on: 30.01.2013
 *      Author: cls
 */

#include "GraphBenchmarkGTest.h"

namespace EnsembleClustering {

GraphBenchmarkGTest::GraphBenchmarkGTest() {
	this->n = 300;
}

GraphBenchmarkGTest::~GraphBenchmarkGTest() {
	// TODO Auto-generated destructor stub
}


// Task: precompute incident weights with different methods


// TEST: use different containers

TEST_F(GraphBenchmarkGTest, useStandard) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	NodeMap<double> incidentWeight(n, 0.0);

	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
	});
	runtime.stop();

	INFO("[DONE] (" << runtime.elapsed().count() << " ms)");

	// test correctness of result
	bool correct = true;
	G.forallNodes([&](node v){
		correct &= (incidentWeight[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}

TEST_F(GraphBenchmarkGTest, useVector) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	std::vector<double> incidentWeight(n+1, 0.0);

	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
	});
	runtime.stop();

	INFO("[DONE] (" << runtime.elapsed().count() << " ms)");

	// test correctness of result
	bool correct = true;
	G.forallNodes([&](node v){
		correct &= (incidentWeight[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}

TEST_F(GraphBenchmarkGTest, useArray) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	double incidentWeight[n + 1];

	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
	});
	runtime.stop();

	INFO("[DONE] (" << runtime.elapsed().count() << " ms)");

	// test correctness of result
	bool correct = true;
	G.forallNodes([&](node v){
		correct &= (incidentWeight[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}


// RESULT: NodeMap, vector and array are about equally fast


// TEST: parallelize with different containers

TEST_F(GraphBenchmarkGTest, useStandard_par) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	NodeMap<double> incidentWeight(n, 0.0);

	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
	}, "parallel");
	runtime.stop();

	INFO("[DONE] (" << runtime.elapsed().count() << " ms)");

	// test correctness of result
	bool correct = true;
	G.forallNodes([&](node v){
		correct &= (incidentWeight[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}

TEST_F(GraphBenchmarkGTest, useVector_par) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	std::vector<double> incidentWeight(n+1, 0.0);

	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
	}, "parallel");
	runtime.stop();

	INFO("[DONE] (" << runtime.elapsed().count() << " ms)");

	// test correctness of result
	bool correct = true;
	G.forallNodes([&](node v){
		correct &= (incidentWeight[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}

TEST_F(GraphBenchmarkGTest, useArray_par) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	double incidentWeight[n + 1];

	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
	}, "parallel");
	runtime.stop();

	INFO("[DONE] (" << runtime.elapsed().count() << " ms)");

	// test correctness of result
	bool correct = true;
	G.forallNodes([&](node v){
		correct &= (incidentWeight[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}



} /* namespace EnsembleClustering */
