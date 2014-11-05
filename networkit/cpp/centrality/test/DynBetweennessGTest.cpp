/*
 * DynBetweennessGTest.cpp
 *
 *  Created on: 05.08.2014
 *      Author: ebergamini, cls
 */

#include "DynBetweennessGTest.h"
#include "../Betweenness.h"
#include "../DynApproxBetweenness.h"
#include "../ApproxBetweenness.h"
#include "../ApproxBetweenness2.h"
#include "../DynBetweenness.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"
#include "../../auxiliary/NumericTools.h"
#include "../../graph/Sampling.h"
#include "../../generators/DorogovtsevMendesGenerator.h"

namespace NetworKit {

TEST_F(DynBetweennessGTest, testDynBetweennessSmallGraph) {
/* Graph:
   0    3   6
	\  / \ /
	2    5
	/  \ / \
   1    4   7
*/
	int n = 8;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);

	DynBetweenness dynbc(G, false);
	Betweenness bc(G);
	dynbc.run();
	bc.run();
	std::vector<double> dynbc_scores = dynbc.scores();
	std::vector<double> bc_scores = bc.scores();

	int i;
	const double tol = 1e-8;
	for(i=0; i<n; i++) {
		EXPECT_NEAR(dynbc_scores[i], bc_scores[i], tol) << "Scores are different";
	}

	// edge insertions
	GraphEvent ev(GraphEvent::EDGE_ADDITION, 0, 6, 1.0);
	G.addEdge(ev.u, ev.v);
	bc.run();
	dynbc.update(ev);

	dynbc_scores = dynbc.scores();
	bc_scores = bc.scores();
	for(i=0; i<n; i++) {
		EXPECT_NEAR(dynbc_scores[i], bc_scores[i], tol) << "Scores are different";
	}

}


TEST_F(DynBetweennessGTest, testWeightedDynBetweennessSmallGraph) {
/* Graph:
0    3   6
	\  / \ /
	2    5
	/  \ / \
1    4   7
*/
	int n = 8;
	Graph G(n, true);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);

	DynBetweenness dynbc(G, true);
//	Graph G1 = Graph(G, false, false);
	Betweenness bc(G);
	dynbc.run();
	bc.run();
	std::vector<double> dynbc_scores = dynbc.scores();
	std::vector<double> bc_scores = bc.scores();

	int i;
	const double tol = 1e-8;
	for(i=0; i<n; i++) {
		EXPECT_NEAR(dynbc_scores[i], bc_scores[i], tol) << "Scores are different";
	}

	// edge insertions
	GraphEvent ev(GraphEvent::EDGE_ADDITION, 0, 6, 1.0);
	G.addEdge(ev.u, ev.v);
//	G1.addEdge(ev.u, ev.v);
	bc.run();
	dynbc.update(ev);
	DEBUG("after");
	dynbc_scores = dynbc.scores();
	bc_scores = bc.scores();
	for(i=0; i<n; i++) {
		EXPECT_NEAR(dynbc_scores[i], bc_scores[i], tol) << "Scores are different";
	}

}

TEST_F(DynBetweennessGTest, testDynApproxBetweennessSmallGraph) {
/* Graph:
0    3   6
	\  / \ /
	2    5
	/  \ / \
1    4   7
*/
	int n = 8;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);

	double epsilon = 0.01; // error
	double delta = 0.1; // confidence
	DynApproxBetweenness dynbc(G, epsilon, delta);
	Betweenness bc(G);
	dynbc.run();
	bc.run();
	std::vector<double> dynbc_scores = dynbc.scores();
	std::vector<double> bc_scores = bc.scores();
	for(int i=0; i<n; i++) {
		DEBUG("Difference ", dynbc_scores[i]-bc_scores[i]/double(n*(n-1)));
	}
	std::vector<GraphEvent> batch;
	batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 0, 6, 1.0));
	G.addEdge(batch[0].u, batch[0].v);
	bc.run();
	dynbc.update(batch);
	dynbc_scores = dynbc.scores();
	bc_scores = bc.scores();
	for(int i=0; i<n; i++) {
		DEBUG("Difference ", dynbc_scores[i]-bc_scores[i]/double(n*(n-1)));
	}

}


TEST_F(DynBetweennessGTest, testDynVsStatic) {
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");
	count n = G.upperNodeIdBound();

	double epsilon = 0.1; // error
	double delta = 0.1; // confidence
	DynApproxBetweenness dynbc(G, epsilon, delta, false);
	ApproxBetweenness bc(G, epsilon, delta);
	dynbc.run();
	bc.run();
	std::vector<double> dynbc_scores = dynbc.scores();
	std::vector<double> bc_scores = bc.scores();
	double err1=0;
	for(count i=0; i<n; i++) {
		double x = dynbc_scores[i]-bc_scores[i];
		if (x > err1)
			err1 = x;
	}
	DEBUG("Before the edge insertion: ");
	std::vector<GraphEvent> batch;
	count nInsertions = 10, i = 0;
	while (i < nInsertions) {
		node v1 = Sampling::randomNode(G);
		node v2 = Sampling::randomNode(G);
		if (v1 != v2 && !G.hasEdge(v1, v2)) {
			G.addEdge(v1, v2);
			batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
			i++;
		}
	}
	bc.run();
	dynbc.update(batch);
	dynbc_scores = dynbc.scores();
	bc_scores = bc.scores();
	err1 = 0;
	for(count i=0; i<n; i++) {
		double x = dynbc_scores[i]-bc_scores[i];
		if (x > err1)
			err1 = x;
	}
	DEBUG("After the edge insertion: ");

}


TEST_F(DynBetweennessGTest, testCorrectnessDynExactBetweenness) {
	METISGraphReader reader;
	DorogovtsevMendesGenerator generator(100);
	Graph G = generator.generate();
	int n = G.upperNodeIdBound();
	DynBetweenness dynbc(G, true);
	Betweenness bc(G);
	dynbc.run();
	bc.run();
	DEBUG("Before the edge insertion: ");
	GraphEvent ev;
	count nInsertions = 1, i = 0;
	while (i < nInsertions) {
		node v1 = Sampling::randomNode(G);
		node v2 = Sampling::randomNode(G);
		if (v1 != v2 && !G.hasEdge(v1, v2)) {
			i++;
			G.addEdge(v1, v2);
			ev = GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0);
			dynbc.update(ev);
			bc.run();
			std::vector<double> dynbc_scores = dynbc.scores();
			std::vector<double> bc_scores = bc.scores();
			int j;
			const double tol = 1e-6;
			for(j=0; j<n; j++) {
				EXPECT_NEAR(dynbc_scores[j], bc_scores[j], tol) << "Scores are different";
			}
		}
	}
}


TEST_F(DynBetweennessGTest, testApproxBetweenness) {
	METISGraphReader reader;
	DorogovtsevMendesGenerator generator(1000);
	Graph G1 = generator.generate();
	Graph G(G1, true, false);
	ApproxBetweenness bc(G, 0.1, 0.1);
	bc.run();
	DEBUG("Number of samples: ", bc.numberOfSamples());
	ApproxBetweenness bc1(G1, 0.1, 0.1);
	bc1.run();
	DEBUG("Number of samples: ", bc1.numberOfSamples());
}


TEST_F(DynBetweennessGTest, testWeightedDynExactBetweenness) {
	METISGraphReader reader;
	DorogovtsevMendesGenerator generator(100);
	Graph G1 = generator.generate();
	Graph G(G1, true, false);
	DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
	int n = G.upperNodeIdBound();
	DynBetweenness dynbc(G, true);
	Betweenness bc(G);
	dynbc.run();
	bc.run();
	DEBUG("Before the edge insertion: ");
	GraphEvent ev;
	count nInsertions = 1, i = 0;
	while (i < nInsertions) {
		DEBUG("Sampling a new edge");
		node v1 = Sampling::randomNode(G);
		node v2 = Sampling::randomNode(G);
		if (v1 != v2 && !G.hasEdge(v1, v2)) {
			i++;
			DEBUG("Adding edge number ", i);
			G.addEdge(v1, v2);
			ev = GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0);
			DEBUG("Running update with dynamic bc");
			dynbc.update(ev);
			DEBUG("Running from scratch with bc");
			bc.run();
			std::vector<double> dynbc_scores = dynbc.scores();
			std::vector<double> bc_scores = bc.scores();
			int j;
			const double tol = 1e-6;
			for(j=0; j<n; j++) {
				EXPECT_NEAR(dynbc_scores[j], bc_scores[j], tol) << "Scores are different";
			}
		}
	}
}

} /* namespace NetworKit */
