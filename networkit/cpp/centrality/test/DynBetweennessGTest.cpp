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

	DynBetweenness dynbc = DynBetweenness(G, false);
	Betweenness bc = Betweenness(G);
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

	DynBetweenness dynbc = DynBetweenness(G, true);
//	Graph G1 = Graph(G, false, false);
	Betweenness bc = Betweenness(G);
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
	DynApproxBetweenness dynbc = DynApproxBetweenness(G, epsilon, delta);
	Betweenness bc = Betweenness(G);
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
	DynApproxBetweenness dynbc = DynApproxBetweenness(G, epsilon, delta, false);
	ApproxBetweenness bc = ApproxBetweenness(G, epsilon, delta);
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

TEST_F(DynBetweennessGTest, timeDynApproxBetweenness) {
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");

	double epsilon = 0.1; // error
	double delta = 0.1; // confidence
	DynApproxBetweenness dynbc = DynApproxBetweenness(G, epsilon, delta, false);
	INFO("initial run");
	dynbc.run();
	INFO("update");
	std::vector<GraphEvent> batch;
	count nInsertions = 100000, i = 0;
	while (i < nInsertions) {
		node v1 = Sampling::randomNode(G);
		node v2 = Sampling::randomNode(G);
		if (v1 != v2 && !G.hasEdge(v1, v2)) {
			G.addEdge(v1, v2);
			batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
			i++;
		}
	}
	dynbc.update(batch);
}

TEST_F(DynBetweennessGTest, timeDynExactBetweenness) {
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");

	DynBetweenness dynbc = DynBetweenness(G);
	INFO("initial run");
	dynbc.run();
	INFO("update");
	GraphEvent ev;
	count nInsertions = 10, i = 0;
	while (i < nInsertions) {
		node v1 = Sampling::randomNode(G);
		node v2 = Sampling::randomNode(G);
		if (v1 != v2 && !G.hasEdge(v1, v2)) {
			G.addEdge(v1, v2);
			ev = GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0);
			dynbc.update(ev);
			i++;
		}
	}
}

TEST_F(DynBetweennessGTest, testCorrectnessDynExactBetweenness) {
	METISGraphReader reader;
	DorogovtsevMendesGenerator generator(1000);
	Graph G = generator.generate();
	int n = G.upperNodeIdBound();
	DynBetweenness dynbc = DynBetweenness(G, true);
	Betweenness bc = Betweenness(G);
	dynbc.run();
	bc.run();
	DEBUG("Before the edge insertion: ");
	GraphEvent ev;
	count nInsertions = 10, i = 0;
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

TEST_F(DynBetweennessGTest, compareAffectedVertices) {
	METISGraphReader reader;
	DorogovtsevMendesGenerator generator(100);
	Graph G = generator.generate();
	DynBetweenness dynbc = DynBetweenness(G, true);
	dynbc.run();
	std::vector<std::vector<edgeweight>> dist1;
	std::vector<std::vector<double>> dep1;
	std::vector<std::vector<bigfloat>> npaths1;
	dist1.resize(G.upperNodeIdBound());
	dep1.resize(G.upperNodeIdBound());
	npaths1.resize(G.upperNodeIdBound());
	G.forNodes([&] (node s){
		dist1[s].resize(G.upperNodeIdBound());
		dep1[s].resize(G.upperNodeIdBound());
		npaths1[s].resize(G.upperNodeIdBound());
		G.forNodes([&] (node t){
			dist1[s][t] = dynbc.distance(s, t);
			dep1[s][t] = dynbc.dependency(s, t);
			npaths1[s][t] = dynbc.nPaths(s, t);
		});
	});
	std::vector<std::vector<edgeweight>> dist2 = dist1;
	std::vector<std::vector<double>> dep2 = dep1;
	auto npaths2 = npaths1;
	DEBUG("Before the edge insertion: ");
	GraphEvent ev;
	count nInsertions = 10, i = 0;
	int totAffectedDep = 0;
	while (i < nInsertions) {
		node v1 = Sampling::randomNode(G);
		node v2 = Sampling::randomNode(G);
		if (v1 != v2 && !G.hasEdge(v1, v2)) {
			i++;
			G.addEdge(v1, v2);
			ev = GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0);
			dynbc.update(ev);
			// compare the old distances, number of shortest paths and dependencies with the new ones
			int diff_dep = 0;
			int diff_dist = 0;
			int diff_dep2 = 0;
			int diff_dist2 = 0;
			G.forNodes([&] (node s){
				G.forNodes([&] (node t){
					if (!Aux::NumericTools::logically_equal(dist1[s][t], dynbc.distance(s,t)) || npaths1[s][t] != dynbc.nPaths(s,t)) {
						diff_dist ++;
					}
					if (!Aux::NumericTools::logically_equal(dist1[s][t], dynbc.distance(s,t)) || (npaths1[s][t] != dynbc.nPaths(s,t)) || !Aux::NumericTools::logically_equal(dep1[s][t], dynbc.dependency(s,t))) {
						diff_dep ++;
					}
					if (!Aux::NumericTools::logically_equal(dist2[s][t], dynbc.distance(s,t)) ||(npaths2[s][t] != dynbc.nPaths(s,t))) {
						diff_dist2 ++;
					}
					if (!Aux::NumericTools::logically_equal(dist2[s][t], dynbc.distance(s,t)) || (npaths2[s][t] != dynbc.nPaths(s,t)) || !Aux::NumericTools::logically_equal(dep2[s][t], dynbc.dependency(s,t))) {
						diff_dep2 ++;
					}
				});
			});
			std::cout<<"Number of vertices affected by the "<<i<<"-th batch"<<std::endl;
			std::cout<<i<<" Diff_dist: "<<diff_dist2<<std::endl;
			std::cout<<i<<" Diff_dep: "<<diff_dep2<<std::endl;
			totAffectedDep += diff_dep2;
			std::cout<<"Total number of vertices affected up to now (counting only once the vertices affected by more than one edge insertion of the batch)"<<std::endl;
			std::cout<<i<<" Diff_dist: "<<diff_dist<<std::endl;
			std::cout<<i<<" Diff_dep: "<<diff_dep<<std::endl;
			dist2.clear();
			dep2.clear();
			npaths2.clear();
			dist2.resize(G.upperNodeIdBound());
			dep2.resize(G.upperNodeIdBound());
			npaths2.resize(G.upperNodeIdBound());
			G.forNodes([&] (node s){
				dist2[s].resize(G.upperNodeIdBound());
				dep2[s].resize(G.upperNodeIdBound());
				npaths2[s].resize(G.upperNodeIdBound());
				G.forNodes([&] (node t){
					dist2[s][t] = dynbc.distance(s, t);
					dep2[s][t] = dynbc.dependency(s, t);
					npaths2[s][t] = dynbc.nPaths(s, t);
				});
			});
		}
	}
	std::cout<<"Sum of vertices whose dependencies have been affected in the batches: "<<totAffectedDep<<std::endl;
}

TEST_F(DynBetweennessGTest, testApproxBetweenness) {
	METISGraphReader reader;
	DorogovtsevMendesGenerator generator(1000);
	Graph G1 = generator.generate();
	Graph G = Graph(G1, true, false);
	ApproxBetweenness bc(G, 0.1, 0.1);
	bc.run();
	DEBUG("Number of samples: ", bc.numberOfSamples());
	ApproxBetweenness bc1(G1, 0.1, 0.1);
	bc1.run();
	DEBUG("Number of samples: ", bc1.numberOfSamples());
}


TEST_F(DynBetweennessGTest, testWeightedDynExactBetweenness) {
	METISGraphReader reader;
	DorogovtsevMendesGenerator generator(1000);
	Graph G1 = generator.generate();
	Graph G = Graph(G1, true, false);
	DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
	int n = G.upperNodeIdBound();
	DynBetweenness dynbc = DynBetweenness(G, true);
	Betweenness bc = Betweenness(G);
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
