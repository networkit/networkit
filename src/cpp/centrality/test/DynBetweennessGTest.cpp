/*
 * DynBetweennessGTest.cpp
 *
 *  Created on: 05.08.2014
 *      Author: ebergamini, cls
 */

#include "DynBetweennessGTest.h"
#include "../Betweenness.h"
#include "../Betweenness2.h"
#include "../DynApproxBetweenness.h"
#include "../ApproxBetweenness.h"
#include "../ApproxBetweenness2.h"
#include "../DynBetweenness.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"
#include "../../graph/Sampling.h"

namespace NetworKit {

TEST_F(DynBetweennessGTest, testDynVsStatic) {
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");
	count n = G.upperNodeIdBound();
	std::cout<<n<<std::endl;

	double epsilon = 0.01; // error
	double delta = 0.1; // confidence
	DynApproxBetweenness dynbc = DynApproxBetweenness(G, epsilon, delta);
	ApproxBetweenness bc = ApproxBetweenness(G, epsilon, delta);
	dynbc.run();
	bc.run();
	std::vector<double> dynbc_scores = dynbc.scores();
	std::vector<double> bc_scores = bc.scores();
	double err1=0;
	for(count i=0; i<n; i++) {
		//std::cout<<dynbc_scores[i]-bc_scores[i]/double(n*(n-1))<<std::endl;
		double x = dynbc_scores[i]-bc_scores[i];
		if (x > err1)
			err1 = x;
	}
	std::cout<<"Before the edge insertion: "<<err1<<std::endl;
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
		//std::cout<<dynbc_scores[i]-bc_scores[i]/double(n*(n-1))<<std::endl;
		double x = dynbc_scores[i]-bc_scores[i];
		if (x > err1)
			err1 = x;
	}
	std::cout<<"After the edge insertion: "<<err1<<std::endl;

}

TEST_F(DynBetweennessGTest, timeDynApproxBetweenness) {
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");

	double epsilon = 0.1; // error
	double delta = 0.1; // confidence
	DynApproxBetweenness dynbc = DynApproxBetweenness(G, epsilon, delta);
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
	count nInsertions = 100, i = 0;
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
	Graph G = reader.read("input/PGPgiantcompo.graph");

	DynBetweenness dynbc = DynBetweenness(G);
	dynbc.run();
	std::cout<<"Before the edge insertion: "<<std::endl;
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


} /* namespace NetworKit */
