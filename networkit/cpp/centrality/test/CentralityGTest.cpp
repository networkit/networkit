/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "CentralityGTest.h"
#include "../Betweenness.h"
#include "../Closeness.h"
#include "../DynApproxBetweenness.h"
#include "../ApproxBetweenness.h"
#include "../EstimateBetweenness.h"
#include "../SpanningEdgeCentrality.h"
#include "../ApproxCloseness.h"
#include "../EigenvectorCentrality.h"
#include "../KatzCentrality.h"
#include "../PageRank.h"
#include "../KPathCentrality.h"
#include "../CoreDecomposition.h"
#include "../LocalClusteringCoefficient.h"
#include "../../io/METISGraphReader.h"
#include "../../io/SNAPGraphReader.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../auxiliary/Log.h"
#include "../../structures/Cover.h"
#include "../PermanenceCentrality.h"
#include "../../structures/Partition.h"
#include "../../auxiliary/Timer.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../generators/DorogovtsevMendesGenerator.h"
#include "../TopCloseness.h"
#include <iostream>
#include <iomanip>



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

TEST_F(CentralityGTest, testKatzCentralityDirected) {
	SNAPGraphReader reader;
	Graph G = reader.read("input/wiki-Vote.txt"); // TODO: replace by smaller graph
	KatzCentrality kc(G);

	DEBUG("start kc run");
	kc.run();
	DEBUG("finish kc");
	std::vector<std::pair<node, double> > kc_ranking = kc.ranking();
	std::vector<double> kc_scores = kc.scores();

	EXPECT_EQ(kc_ranking[0].first, 699u);
}

TEST_F(CentralityGTest, testPageRankDirected) {
	SNAPGraphReader reader;
	Graph G = reader.read("input/wiki-Vote.txt"); // TODO: replace by smaller graph
	PageRank pr(G);

	DEBUG("start pr run");
	pr.run();
	DEBUG("finish pr");
	std::vector<std::pair<node, double> > pr_ranking = pr.ranking();

	const double tol = 1e-3;
	EXPECT_EQ(pr_ranking[0].first, 699u);
	EXPECT_NEAR(pr_ranking[0].second, 0.00432, tol);
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



TEST_F(CentralityGTest, testEstimateBetweenness) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");

	EstimateBetweenness abc2(G, 100);
	abc2.run();

	DEBUG("approximated betweenness scores: ", abc2.scores());
}

TEST_F(CentralityGTest, testApproxClosenessCentralityOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");

	ApproxCloseness acc(G, 100);
	acc.run();

	DEBUG("approximated closeness scores: ", acc.scores());
}

TEST_F(CentralityGTest, testApproxClosenessCentralityOnToyGraph) {
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

    ApproxCloseness acc(G, 6, 0.1, false);
    acc.run();
    std::vector<double> cc = acc.scores();

	double maximum = acc.maximum();

    const double tol = 0.2;
    EXPECT_NEAR(0.1, cc[0], tol);
    EXPECT_NEAR(0.1, cc[1], tol);
    EXPECT_NEAR(0.166667, cc[2], tol);
    EXPECT_NEAR(0.125, cc[3], tol);
    EXPECT_NEAR(0.125, cc[4], tol);
    EXPECT_NEAR(0.1, cc[5], tol);
	EXPECT_NEAR(0.2, maximum, tol);

	ApproxCloseness acc2(G, 4, 0.1, true);
	acc2.run();
	std::vector<double> cc2 = acc2.scores();

	double maximum2 = acc2.maximum();

	EXPECT_NEAR(0.5, cc2[0], tol);
	EXPECT_NEAR(0.5, cc2[1], tol);
	EXPECT_NEAR(0.833335, cc2[2], tol);
	EXPECT_NEAR(0.625, cc2[3], tol);
	EXPECT_NEAR(0.625, cc2[4], tol);
	EXPECT_NEAR(0.5, cc2[5], tol);
	EXPECT_NEAR(0.2, maximum2, tol);
}

TEST_F(CentralityGTest, testEdgeBetweennessCentrality) {
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
	G.indexEdges();

	Betweenness centrality(G,false,true);
	centrality.run();
	std::vector<double> bc = centrality.edgeScores();

	const double tol = 1e-3;
	EXPECT_NEAR(10.0, bc[0], tol);
	EXPECT_NEAR(10.0, bc[1], tol);
	EXPECT_NEAR(10.0, bc[2], tol);
	EXPECT_NEAR(10.0, bc[3], tol);
	EXPECT_NEAR(6.0, bc[4], tol);
	EXPECT_NEAR(6.0, bc[5], tol);
}


TEST_F(CentralityGTest, tryEdgeBetweennessCentrality) {
    auto path = "input/PGPgiantcompo.graph";
    METISGraphReader reader;
    Graph G = reader.read(path);
    G.indexEdges();

	Betweenness centrality(G,false,true);
	centrality.run();
	std::vector<double> bc = centrality.edgeScores();

}



TEST_F(CentralityGTest, testClosenessCentrality) {
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

    Closeness centrality(G,false);
    centrality.run();
    std::vector<double> bc = centrality.scores();

		double maximum = centrality.maximum();

    const double tol = 1e-3;
    EXPECT_NEAR(0.1, bc[0], tol);
    EXPECT_NEAR(0.1, bc[1], tol);
    EXPECT_NEAR(0.166667, bc[2], tol);
    EXPECT_NEAR(0.125, bc[3], tol);
    EXPECT_NEAR(0.125, bc[4], tol);
    EXPECT_NEAR(0.1, bc[5], tol);
		EXPECT_NEAR(0.2, maximum, tol);
}


TEST_F(CentralityGTest, testKPathCentrality) {
    METISGraphReader reader;
    Graph G = reader.read("input/power.graph");

    KPathCentrality centrality(G);
    centrality.run();
}


TEST_F(CentralityGTest, testCoreDecompositionSimple) {
	count n = 3;
	Graph G(n);
	G.addEdge(0,1);

	CoreDecomposition coreDec(G);
	coreDec.run();
	std::vector<double> coreness = coreDec.scores();

	EXPECT_EQ(1u, coreness[0]) << "expected coreness";
	EXPECT_EQ(1u, coreness[1]) << "expected coreness";
	EXPECT_EQ(0u, coreness[2]) << "expected coreness";
}

TEST_F(CentralityGTest, testCoreDecomposition) {
	count n = 16;
	Graph G(n);

// 	// create graph used in Baur et al. and network analysis lecture
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 5);
	G.addEdge(5, 7);
	G.addEdge(6, 7);

	G.addEdge(6, 8);
	G.addEdge(6, 9);
	G.addEdge(6, 11);
	G.addEdge(7, 12);
	G.addEdge(8, 9);

	G.addEdge(8, 10);
	G.addEdge(8, 11);
	G.addEdge(8, 13);
	G.addEdge(9, 10);
	G.addEdge(9, 11);

	G.addEdge(9, 13);
	G.addEdge(10, 11);
	G.addEdge(10, 13);
	G.addEdge(10, 14);
	G.addEdge(11, 13);

	G.addEdge(11, 14);
	G.addEdge(12, 15);
	G.addEdge(13, 14);
	G.addEdge(14, 15);

	EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
	EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

	// compute core decomposition
	CoreDecomposition coreDec(G);
	coreDec.run();
	std::vector<double> coreness = coreDec.scores();
	// init cores
	// init shells
	Cover cores = coreDec.getCover();
	Partition shells = coreDec.getPartition();

	EXPECT_EQ(0u, coreness[0]) << "expected coreness";
	EXPECT_EQ(0u, coreness[1]) << "expected coreness";
	EXPECT_EQ(1u, coreness[2]) << "expected coreness";
	EXPECT_EQ(1u, coreness[3]) << "expected coreness";
	EXPECT_EQ(1u, coreness[4]) << "expected coreness";
	EXPECT_EQ(1u, coreness[5]) << "expected coreness";
	EXPECT_EQ(3u, coreness[6]) << "expected coreness";
	EXPECT_EQ(2u, coreness[7]) << "expected coreness";
	EXPECT_EQ(4u, coreness[8]) << "expected coreness";
	EXPECT_EQ(4u, coreness[9]) << "expected coreness";
	EXPECT_EQ(4u, coreness[10]) << "expected coreness";
	EXPECT_EQ(4u, coreness[11]) << "expected coreness";
	EXPECT_EQ(2u, coreness[12]) << "expected coreness";
	EXPECT_EQ(4u, coreness[13]) << "expected coreness";
	EXPECT_EQ(3u, coreness[14]) << "expected coreness";
	EXPECT_EQ(2u, coreness[15]) << "expected coreness";

	// for (index e = 0; e < n; e++) {
	// 	EXPECT_EQ(cores.contains(e), true);
	// 	EXPECT_EQ(shells.contains(e), true);
	// }
	// EXPECT_EQ(cores.get, coreness[15]) << "expected coreness";

	// test throw runtime error for self-loop in graph
	Graph H(2);
	H.addEdge(0, 1);
	H.addEdge(1, 1);
	EXPECT_ANY_THROW(CoreDecomposition CoreDec(H));
}

TEST_F(CentralityGTest, benchCoreDecompositionSnapGraphs) {
	SNAPGraphReader reader;
	std::vector<std::string> filenames = {"soc-LiveJournal1.edgelist-t0.graph", "cit-Patents.txt", "com-orkut.ungraph.txt", "web-BerkStan.edgelist-t0.graph"};

	for (auto f: filenames) {
		std::string filename("/home/i11/staudt/Graphs/Static/SNAP/" + f);
		DEBUG("about to read file ", filename);
		Graph G = reader.read(filename);
		G.removeSelfLoops();
		CoreDecomposition coreDec(G, false);
		Aux::Timer timer;
		timer.start();
		coreDec.run();
		timer.stop();
		INFO("Time for ParK of ", filename, ": ", timer.elapsedTag());

		CoreDecomposition coreDec2(G, true);
		timer.start();
		coreDec2.run();
		timer.stop();
		INFO("Time for bucket queue based k-core decomposition of ", filename, ": ", timer.elapsedTag());

		G.forNodes([&](node u) {
			EXPECT_EQ(coreDec.score(u), coreDec2.score(u));
		});
	}
}

TEST_F(CentralityGTest, benchCoreDecompositionDimacsGraphs) {
  METISGraphReader reader;
  std::vector<std::string> filenames = {"eu-2005", "coPapersDBLP", "uk-2002", "europe-osm", "cage15"};

  for (auto f: filenames) {
    std::string filename("/home/i11/staudt/Graphs/Static/DIMACS/Large/" + f + ".metis.graph");
    DEBUG("about to read file ", filename);
    Graph G = reader.read(filename);
    G.removeSelfLoops();
    CoreDecomposition coreDec(G, false);
    Aux::Timer timer;
    timer.start();
    coreDec.run();
    timer.stop();
    INFO("Time for ParK of ", filename, ": ", timer.elapsedTag());

    CoreDecomposition coreDec2(G, true);
    timer.start();
    coreDec2.run();
    timer.stop();
    INFO("Time for bucket queue based k-core decomposition of ", filename, ": ", timer.elapsedTag());

    G.forNodes([&](node u) {
	EXPECT_EQ(coreDec.score(u), coreDec2.score(u));
      });
  }
}

TEST_F(CentralityGTest, benchCoreDecompositionLocal) {
  METISGraphReader reader;
  std::vector<std::string> filenames = {"coPapersCiteseer", "in-2004", "coAuthorsDBLP", "audikw1"};

  for (auto f: filenames) {
    std::string filename("input/" + f + ".graph");
    DEBUG("about to read file ", filename);
    Graph G = reader.read(filename);
    G.removeSelfLoops();
    CoreDecomposition coreDec(G, false);
    Aux::Timer timer;
    timer.start();
    coreDec.run();
    timer.stop();
    INFO("Time for ParK of ", filename, ": ", timer.elapsedTag());

    CoreDecomposition coreDec2(G, true);
    timer.start();
    coreDec2.run();
    timer.stop();
    INFO("Time for bucket queue based k-core decomposition of ", filename, ": ", timer.elapsedTag());

    G.forNodes([&](node u) {
	EXPECT_EQ(coreDec.score(u), coreDec2.score(u));
      });
  }
}

TEST_F(CentralityGTest, testCoreDecompositionDirected) {
	count n = 16;
	Graph G(n,false,true);

// 	// create graph used in Baur et al. and network analysis lecture
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 5);
	G.addEdge(5, 7);
	G.addEdge(6, 7);

	G.addEdge(6, 8);
	G.addEdge(6, 9);
	G.addEdge(6, 11);
	G.addEdge(7, 12);
	G.addEdge(8, 9);

	G.addEdge(8, 10);
	G.addEdge(8, 11);
	G.addEdge(8, 13);
	G.addEdge(9, 10);
	G.addEdge(9, 11);

	G.addEdge(9, 13);
	G.addEdge(10, 11);
	G.addEdge(10, 13);
	G.addEdge(10, 14);
	G.addEdge(11, 13);

	G.addEdge(11, 14);
	G.addEdge(12, 15);
	G.addEdge(13, 14);
	G.addEdge(14, 15);

	EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
	EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

	// compute core decomposition
	CoreDecomposition coreDec(G);
	coreDec.run();
	std::vector<double> coreness = coreDec.scores();

	EXPECT_EQ(0u, coreness[0]) << "expected coreness";
	EXPECT_EQ(0u, coreness[1]) << "expected coreness";
	EXPECT_EQ(1u, coreness[2]) << "expected coreness";
	EXPECT_EQ(1u, coreness[3]) << "expected coreness";
	EXPECT_EQ(1u, coreness[4]) << "expected coreness";
	EXPECT_EQ(1u, coreness[5]) << "expected coreness";
	EXPECT_EQ(3u, coreness[6]) << "expected coreness";
	EXPECT_EQ(2u, coreness[7]) << "expected coreness";
	EXPECT_EQ(4u, coreness[8]) << "expected coreness";
	EXPECT_EQ(4u, coreness[9]) << "expected coreness";
	EXPECT_EQ(4u, coreness[10]) << "expected coreness";
	EXPECT_EQ(4u, coreness[11]) << "expected coreness";
	EXPECT_EQ(2u, coreness[12]) << "expected coreness";
	EXPECT_EQ(4u, coreness[13]) << "expected coreness";
	EXPECT_EQ(3u, coreness[14]) << "expected coreness";
	EXPECT_EQ(2u, coreness[15]) << "expected coreness";
}

TEST_F(CentralityGTest, testLocalClusteringCoefficientUndirected) {
	count n = 16;
	Graph G(n,false,false);

// 	// create graph used in Baur et al. and network analysis lecture
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 5);
	G.addEdge(5, 7);
	G.addEdge(6, 7);

	G.addEdge(6, 8);
	G.addEdge(6, 9);
	G.addEdge(6, 11);
	G.addEdge(7, 12);
	G.addEdge(8, 9);

	G.addEdge(8, 10);
	G.addEdge(8, 11);
	G.addEdge(8, 13);
	G.addEdge(9, 10);
	G.addEdge(9, 11);

	G.addEdge(9, 13);
	G.addEdge(10, 11);
	G.addEdge(10, 13);
	G.addEdge(10, 14);
	G.addEdge(11, 13);

	G.addEdge(11, 14);
	G.addEdge(12, 15);
	G.addEdge(13, 14);
	G.addEdge(14, 15);

	EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
	EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

	// compute core decomposition
	LocalClusteringCoefficient lcc(G);
	lcc.run();
	std::vector<double> lccScores = lcc.scores();
	std::vector<double> reference = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.8, 0.8, 0.8, 0.6666666666666666, 0.0, 0.8, 0.5, 0.0};

 	EXPECT_EQ(reference,lccScores);

	LocalClusteringCoefficient lccTurbo(G, true);
	lccTurbo.run();
	EXPECT_EQ(reference, lccTurbo.scores());

	// test throw runtime error for self-loop in graph
	Graph H(2);
	H.addEdge(0, 1);
	H.addEdge(1, 1);
	EXPECT_ANY_THROW(LocalClusteringCoefficient lcc(H));
}

TEST_F(CentralityGTest, testLocalClusteringCoefficientUndirected2) {
	Graph G(6,false,false);
	G.addEdge(1, 0);
	G.addEdge(2, 0);
	G.addEdge(2, 1);
	G.addEdge(3, 2);
	G.addEdge(3, 0);
	G.addEdge(3, 1);
	G.addEdge(4, 2);
	G.addEdge(4, 0);
	G.addEdge(5, 3);
	G.addEdge(5, 4);
	G.addEdge(5, 1);
	LocalClusteringCoefficient lcc(G);
	lcc.run();
	std::vector<double> lccScores = lcc.scores();
	std::vector<double> reference = {0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.3333333333333333, 0.3333333333333333};

 	EXPECT_EQ(reference,lccScores);
 }

TEST_F(CentralityGTest, testSimplePermanence) {
	Graph G(15, false, false);
	G.addEdge(0, 1);
	G.addEdge(1, 2);
	G.addEdge(2, 0);
	G.addEdge(2, 3);
	node v = 4;
	node u = 5;
	G.addEdge(v, 0);
	G.addEdge(v, 1);
	G.addEdge(v, 2);
	G.addEdge(u, 3);
	G.addEdge(u, 2);
	G.addEdge(u, 0);
	G.addEdge(6, 7);
	G.addEdge(7, 8);
	G.addEdge(u, 6);
	G.addEdge(u, 7);
	G.addEdge(u, 8);
	G.addEdge(v, 6);
	G.addEdge(v, 7);
	G.addEdge(9, 10);
	G.addEdge(10, 11);
	G.addEdge(u, 9);
	G.addEdge(v, 10);
	G.addEdge(v, 11);
	G.addEdge(12, 13);
	G.addEdge(13, 14);
	G.addEdge(12, 14);
	G.addEdge(v, 12);
	G.addEdge(v, 14);

	Partition P(G.upperNodeIdBound());
	P.setUpperBound(4);
	P[0] = 0;
	P[1] = 0;
	P[2] = 0;
	P[3] = 0;
	P[v] = 0;
	P[u] = 0;
	P[6] = 1;
	P[7] = 1;
	P[8] = 1;
	P[9] = 2;
	P[10] = 2;
	P[11] = 2;
	P[12] = 3;
	P[13] = 3;
	P[14] = 3;

	ASSERT_EQ(9u, G.degree(v));
	ASSERT_EQ(7u, G.degree(u));

	PermanenceCentrality perm(G, P);
	perm.run();
	EXPECT_DOUBLE_EQ(2.0/3.0, perm.getIntraClustering(u));
	EXPECT_DOUBLE_EQ(1, perm.getIntraClustering(v));

	EXPECT_NEAR(-0.19048, perm.getPermanence(u), 0.0005);
	EXPECT_NEAR(0.167, perm.getPermanence(v), 0.0005);
}

TEST_F(CentralityGTest, testTopClosenessDirected) {
    Aux::Random::setSeed(42, true);
    count size = 400;
    count k = 10;
    Graph G1 = DorogovtsevMendesGenerator(size).generate();
    Graph G(G1.upperNodeIdBound(), false, true);
    G1.forEdges([&](node u, node v){
        G.addEdge(u, v);
        G.addEdge(v, u);
    });
    INFO("Number of nodes: ", G.upperNodeIdBound(), ", number of edges: ", G.numberOfEdges());
    Closeness cc(G1, true);
    cc.run();
    TopCloseness topcc(G, k, true, true);
    topcc.run();
    const edgeweight tol = 1e-7;
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc.topkScoresList()[i], tol);
    }
    TopCloseness topcc2(G, k, true, false);
    topcc2.run();
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc2.topkScoresList()[i], tol);
    }
}

TEST_F(CentralityGTest, testTopClosenessUndirected) {
    Aux::Random::setSeed(42, true);
    count size = 400;
    count k = 10;
    Graph G1 = DorogovtsevMendesGenerator(size).generate();
    Graph G(G1.upperNodeIdBound(), false, false);
    G1.forEdges([&](node u, node v){
        G.addEdge(u, v);
        G.addEdge(v, u);
    });
    INFO("Number of nodes: ", G.upperNodeIdBound(), ", number of edges: ", G.numberOfEdges());
    Closeness cc(G1, true);
    cc.run();
    TopCloseness topcc(G, k, true, true);
    topcc.run();
    const edgeweight tol = 1e-7;
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc.topkScoresList()[i], tol);
    }
    TopCloseness topcc2(G, k, true, false);
    topcc2.run();
    for (count i = 0; i < k; i++) {
        EXPECT_NEAR(cc.ranking()[i].second, topcc2.topkScoresList()[i], tol);
    }
}

} /* namespace NetworKit */
