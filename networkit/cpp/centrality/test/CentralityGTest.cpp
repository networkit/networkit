/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "CentralityGTest.h"
#include "../Betweenness.h"
#include "../DynApproxBetweenness.h"
#include "../ApproxBetweenness.h"
#include "../ApproxBetweenness2.h"
#include "../EigenvectorCentrality.h"
#include "../PageRank.h"
#include "../DynBetweenness.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"

namespace NetworKit {

TEST_F(CentralityGTest, testBetweennessCentrality) {
 /* Graph:
    0    3
     \  / \
      2    5
     /  \ /
    1    4
 */
	count n = 6;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);

	Betweenness centrality(G);
	centrality.run();
	std::vector<double> bc = centrality.scores();

	const double tol = 1e-3;
	EXPECT_NEAR(0.0, bc[0], tol);
	EXPECT_NEAR(0.0, bc[1], tol);
	EXPECT_NEAR(15.0, bc[2], tol);
	EXPECT_NEAR(3.0, bc[3], tol);
	EXPECT_NEAR(3.0, bc[4], tol);
	EXPECT_NEAR(1.0, bc[5], tol);
}


TEST_F(CentralityGTest, testBetweenness2Centrality) {
/* Graph:
	0    3
	\  / \
	2    5
	/  \ /
	1    4
*/
	count n = 6;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);

	Betweenness centrality(G);
	centrality.run();
	std::vector<double> bc = centrality.scores();

	const double tol = 1e-3;
	EXPECT_NEAR(0.0, bc[0], tol);
	EXPECT_NEAR(0.0, bc[1], tol);
	EXPECT_NEAR(15.0, bc[2], tol);
	EXPECT_NEAR(3.0, bc[3], tol);
	EXPECT_NEAR(3.0, bc[4], tol);
	EXPECT_NEAR(1.0, bc[5], tol);
}


TEST_F(CentralityGTest, testApproxBetweennessSmallGraph) {
 /* Graph:
    0    3
     \  / \
      2    5
     /  \ /
    1    4
 */
	count n = 6;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);

	double epsilon = 0.1; // error
	double delta = 0.1; // confidence
	ApproxBetweenness centrality(G, epsilon, delta, 0);
	centrality.run();
	std::vector<double> bc = centrality.scores();

	DEBUG("scores: ", bc);
}

TEST_F(CentralityGTest, tryApproxBetweennessOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/ns894786.mps.gz.variable.graph");

	double epsilon = 0.01; // error
	double delta = 0.1; // confidence
	ApproxBetweenness centrality(G, epsilon, delta);
	centrality.run();
	std::vector<double> bc = centrality.scores();

	DEBUG("scores: ", bc);
}


TEST_F(CentralityGTest, testBetweennessCentralityWeighted) {
 /* Graph:
    0    3   6
     \  / \ /
      2 -- 5
     /  \ / \
    1    4   7

    Edges in the upper row have weight 3,
    the edge in the middle row has weight 1.5,
    edges in the lower row have weight 2.
 */
	count n = 8;
	Graph G(n, true);

	G.addEdge(0, 2, 3);
	G.addEdge(1, 2, 2);
	G.addEdge(2, 3, 3);
	G.addEdge(2, 4, 2);
	G.addEdge(2, 5, 1.5);
	G.addEdge(3, 5, 3);
	G.addEdge(4, 5, 2);
	G.addEdge(5, 6, 3);
	G.addEdge(5, 7, 2);

	Betweenness centrality(G);
	centrality.run();
	std::vector<double> bc = centrality.scores();

	const double tol = 1e-3;
	EXPECT_NEAR(0.0, bc[0], tol);
	EXPECT_NEAR(0.0, bc[1], tol);
	EXPECT_NEAR(23.0, bc[2], tol);
	EXPECT_NEAR(0.0, bc[3], tol);
	EXPECT_NEAR(0.0, bc[4], tol);
	EXPECT_NEAR(23.0, bc[5], tol);
	EXPECT_NEAR(0.0, bc[6], tol);
	EXPECT_NEAR(0.0, bc[7], tol);
}

TEST_F(CentralityGTest, testEigenvectorCentrality) {
 /* Graph:
    0    3   6
     \  / \ /
      2 -- 5
     /  \ / \
    1    4   7

    Edges in the upper row have weight 3,
    the edge in the middle row has weight 1.5,
    edges in the lower row have weight 2.
 */
	count n = 8;
	Graph G(n, true);

	G.addEdge(0, 2, 3);
	G.addEdge(1, 2, 2);
	G.addEdge(2, 3, 3);
	G.addEdge(2, 4, 2);
	G.addEdge(2, 5, 1.5);
	G.addEdge(3, 5, 3);
	G.addEdge(4, 5, 2);
	G.addEdge(5, 6, 3);
	G.addEdge(5, 7, 2);

	EigenvectorCentrality centrality(G);
	centrality.run();
	std::vector<double> cen = centrality.scores();

	// computed with Matlab
	const double tol = 1e-4;
	EXPECT_NEAR(0.2254, fabs(cen[0]), tol);
	EXPECT_NEAR(0.1503, fabs(cen[1]), tol);
	EXPECT_NEAR(0.5290, fabs(cen[2]), tol);
	EXPECT_NEAR(0.4508, fabs(cen[3]), tol);
	EXPECT_NEAR(0.3006, fabs(cen[4]), tol);
	EXPECT_NEAR(0.5290, fabs(cen[5]), tol);
	EXPECT_NEAR(0.2254, fabs(cen[6]), tol);
	EXPECT_NEAR(0.1503, fabs(cen[7]), tol);
}

TEST_F(CentralityGTest, testPageRankCentrality) {
 /* Graph:
    0    3   6
     \  / \ /
      2 -- 5
     /  \ / \
    1    4   7

    Edges in the upper row have weight 3,
    the edge in the middle row has weight 1.5,
    edges in the lower row have weight 2.
 */
	count n = 8;
	Graph G(n, true);

	G.addEdge(0, 2, 3);
	G.addEdge(1, 2, 2);
	G.addEdge(2, 3, 3);
	G.addEdge(2, 4, 2);
	G.addEdge(2, 5, 1.5);
	G.addEdge(3, 5, 3);
	G.addEdge(4, 5, 2);
	G.addEdge(5, 6, 3);
	G.addEdge(5, 7, 2);

	double damp = 0.85;
	PageRank centrality(G, damp);
	centrality.run();
	std::vector<double> cen = centrality.scores();

	// compare to Matlab results
	const double tol = 1e-4;
	EXPECT_NEAR(0.0753, fabs(cen[0]), tol);
	EXPECT_NEAR(0.0565, fabs(cen[1]), tol);
	EXPECT_NEAR(0.2552, fabs(cen[2]), tol);
	EXPECT_NEAR(0.1319, fabs(cen[3]), tol);
	EXPECT_NEAR(0.0942, fabs(cen[4]), tol);
	EXPECT_NEAR(0.2552, fabs(cen[5]), tol);
	EXPECT_NEAR(0.0753, fabs(cen[6]), tol);
	EXPECT_NEAR(0.0565, fabs(cen[7]), tol);
}

TEST_F(CentralityGTest, benchSequentialBetweennessCentralityOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");
	Betweenness bc(G);
	bc.run();
	std::vector<std::pair<node, double> > ranking = bc.ranking();
	INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchParallelBetweennessCentralityOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");
	Betweenness bc(G);
	bc.run();
	std::vector<std::pair<node, double> > ranking = bc.ranking();
	INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchEigenvectorCentralityOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");
	EigenvectorCentrality cen(G);
	cen.run();
	std::vector<std::pair<node, double> > ranking = cen.ranking();
	INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchPageRankCentralityOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");
	double damp = 0.85;
	PageRank cen(G, damp);
	cen.run();
	std::vector<std::pair<node, double> > ranking = cen.ranking();
	INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}



TEST_F(CentralityGTest, testApproxBetweenness2) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");

	ApproxBetweenness2 abc2(G, 100);
	abc2.run();

	DEBUG("approximated betweenness scores: ", abc2.scores());

}


} /* namespace NetworKit */
